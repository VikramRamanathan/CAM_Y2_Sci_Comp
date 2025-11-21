import numpy as np
from enum import Enum

# Simulates 2D collisions, because I enjoy pain, I suppose
# This is effectively a thermal expansion simulator, with the right hand side wall of the container
# being a piston with variable mass, with a variable external atmospheric pressure, and the other walls
# being temperature reserviors of some thermal conductivity and temperature.


def compute_cell_from_position(position, grid_cell_height, grid_cell_width):
    return int(position[0] // grid_cell_width), int(position[1] // grid_cell_height)

class Collided(Enum):
    COLLISION = 1
    NO_COLLISION = 2


# Not constant, allows for future insertion
num_grid_rows = 8
num_grid_cols = 8
arraygrid = [[[] for i in range(num_grid_rows)] for j in range(num_grid_cols)] # The arraygrid is indexed x, y

# The corners of the grid
grid_bounds = [0, 0, 7, 7]

# These are consts, that together with the number of rows define the size of the space
grid_cell_height = 64
grid_cell_width = 64

# Using this class is going to add quite a lot of overhead
# But it's... fine, probably
class Particle:
    position: list
    velocity: list
    mass: float
    radius: float
    cell_coords: list
    
    def move_particle(self, new_pos, grid_cell_height, grid_cell_width, num_grid_cols, num_grid_rows, arraygrid):
        new_cell_coords = compute_cell_from_position(new_pos, grid_cell_height, grid_cell_width)
        self.position = new_pos
        
        if self.cell_coords != new_cell_coords:
            curr_cell = arraygrid[self.cell_coords[0] - grid_bounds[0]][self.cell_coords[1] - grid_bounds[1]]
            
            # TODO test this thing
            if new_cell_coords[0] < grid_bounds[0]:
                
                num_new_cols = grid_bounds[0] - new_cell_coords[0]
                num_grid_cols = num_grid_cols + num_new_cols
                new_cols = [[[] for i in range(num_grid_rows)] for column in range(num_new_cols)]
                arraygrid = new_cols + arraygrid
                
                grid_bounds[0] = grid_bounds[0] - num_new_cols
                
            elif new_cell_coords[0] > grid_bounds[2]:
                
                num_new_cols = new_cell_coords[0] - grid_bounds[2]
                num_grid_cols = num_grid_cols + num_new_cols
                new_cols = [[[] for i in range(num_grid_rows)] for column in range(num_new_cols)]
                arraygrid = arraygrid + new_cols
                
                grid_bounds[2] = grid_bounds[2] + num_new_cols
            
            if new_cell_coords[1] < grid_bounds[1]:
                num_new_rows = grid_bounds[1] - new_cell_coords[1]
                num_grid_rows = num_grid_rows + num_new_rows
                arraygrid = [([[]] * num_new_rows) + column for column in arraygrid]
                
                grid_bounds[1] = grid_bounds[1] + num_new_rows
                
            elif new_cell_coords[1] > grid_bounds[3]:
                num_new_rows = new_cell_coords[1] - grid_bounds[3]
                num_grid_rows = num_grid_rows + num_new_rows
                arraygrid = [column + ([[]] * num_new_rows) for column in arraygrid]
                
                grid_bounds[3] = grid_bounds[3] + num_new_rows
            
            #print("Arraygrid len: ", len(arraygrid))
            
            for i in range(len(curr_cell)):
                if curr_cell[i] is self:
                    arraygrid[new_cell_coords[0] - grid_bounds[0]][new_cell_coords[1] - grid_bounds[1]].append(self)
                    self.cell_coords = new_cell_coords
                    curr_cell.pop(i)
                    return num_grid_cols, num_grid_rows, arraygrid
            raise RuntimeError("Particle could not be found in the arraygrid")
        
        return num_grid_cols, num_grid_rows, arraygrid
        
    
    # Does not automatically insert itself into the arraygrid
    def __init__(self, position, velocity, mass, radius, grid_cell_height, grid_cell_width):
        cell_coords = compute_cell_from_position(position, grid_cell_height, grid_cell_width)
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.radius = radius
        self.cell_coords = cell_coords

particles = []
num_particles = 25
avg_init_velocities = 12
mass = 1
radius = 20
for i in range(num_particles):
    pos = [np.random.uniform() * grid_cell_width * num_grid_cols, np.random.uniform() * grid_cell_height * num_grid_rows]
    particle = Particle(pos, np.random.uniform(-1, 1, 2) * 2 * avg_init_velocities, mass, radius, grid_cell_height, grid_cell_width)
    #arraygrid[particle.cell_coords[0] - grid_bounds[0]][particle.cell_coords[1] - grid_bounds[1]].append(particle)
    
    # Lists store things by reference anyway, so it's fine
    particles.append(particle)
    
time = 0
timestep_size = 0.02

# List of 4-tuples, where each tuple defines two pairs of (x, y) coordinates describing the points which the wall lies on
walls = [[0, 0, 0, num_grid_rows * grid_cell_height], [0, 0, 99999, 0], [0, num_grid_rows * grid_cell_height, 99999, num_grid_rows * grid_cell_height]]

# 7 tuple, with the 5th number being mass and the 6th and 7th being velocity components. Rotation is not allowed
moveable_walls = [[num_grid_cols * grid_cell_width, 0, num_grid_cols * grid_cell_width, num_grid_rows * grid_cell_height, 500, 0, 0]]

external_pressure = 100

def clamp(n, smallest, largest):
    return max(smallest, min(n, largest))

def get_cells_to_check(x_coord, y_coord, arraygrid, min_x, min_y, max_x, max_y):
    
    upper_bound_x = clamp(x_coord + 1, min_x, max_x) - grid_bounds[0]
    lower_bound_x = clamp(x_coord - 1, min_x, max_x) - grid_bounds[0]
    upper_bound_y = clamp(y_coord + 1, min_y, max_y) - grid_bounds[1]
    lower_bound_y = clamp(y_coord - 1, min_y, max_y) - grid_bounds[1]
    
    return [[x, y] for x in range(lower_bound_x, upper_bound_x + 1) for y in range(lower_bound_y, upper_bound_y + 1)]
        
    
# Test using numpy arrays for the mass computation of position, velocity and radius
# Overhead of initialisation is difficult to guess at
def check_for_double_particle_collisions(particle, particles_to_check, num_grid_cols, num_grid_rows, arraygrid, timestep_size):
    collisions = []
    
    # I make the very bold assumption that each particle will only collide with at most one other particle in one timestep
    for particle_2 in particles_to_check:
        x_dist = particle_2.position[0] - particle.position[0]
        y_dist = particle_2.position[1] - particle.position[1]
        
        x_movement = (particle_2.velocity[0] - particle.velocity[0]) * timestep_size
        y_movement = (particle_2.velocity[1] - particle.velocity[1]) * timestep_size
        
        
        closest_approach_parameter = (x_dist * x_movement + y_dist * y_movement)
        
        # If the dot product is positive, i.e. the second particle is moving away from our particle, then it cannot hit us
        if closest_approach_parameter >= 0:
            continue
        
        # Cap the movement. This parameter corresponds directly to time
        movement_mag_squared = (x_movement ** 2 + y_movement ** 2)
        inflated_radius_squared = (particle_2.radius + particle.radius) ** 2
        
        # Skip if impossible to reach. Tragically uses a square root to evaluate this, but avoiding that's a micro optimisation at this point
        if movement_mag_squared + inflated_radius_squared + 2 * (movement_mag_squared * inflated_radius_squared)**0.5 < x_dist ** 2 + y_dist ** 2:
            continue
            
        closest_approach_parameter = min(-closest_approach_parameter / movement_mag_squared, 1) # You need to multiply by the normal twice, hence mag squared
        
        closest_approach_x = x_dist + x_movement * closest_approach_parameter
        closest_approach_y = y_dist + y_movement * closest_approach_parameter
        
        closest_approach_dist_squared = closest_approach_x**2 + closest_approach_y**2
        
        if closest_approach_dist_squared <= inflated_radius_squared:
            # They have collided, keep going to find the first collision
            collisions.append((closest_approach_parameter, particle_2, x_dist, y_dist, x_movement, y_movement, closest_approach_x, closest_approach_y,
            closest_approach_dist_squared, inflated_radius_squared))
    
    # No collisions, we are done
    if len(collisions) == 0:
        return arraygrid, Collided.NO_COLLISION
    
    # Find the first collision that will take place
    min_param = 1
    curr_selected = 0
    for i in range(len(collisions)):
        if collisions[i][0] <= min_param:
            curr_selected = i
            min_param = collisions[i][0]
    
    closest_approach_parameter, particle_2, x_dist, y_dist, x_movement, y_movement, closest_approach_x, closest_approach_y, closest_approach_dist_squared, inflated_radius_squared = collisions[i]
    
    # The intersection point is given by a polynomial
    a = x_movement * x_movement + y_movement * y_movement
    b = - x_dist * x_movement - y_dist * y_movement
    c = x_dist * x_dist + y_dist * y_dist - inflated_radius_squared
    
    sqrt_disc = (b**2 - a * c) ** 0.5
    point_1_param = abs((-b - sqrt_disc) / a)
    point_2_param = abs((-b + sqrt_disc) / a)
    
    if point_1_param < point_2_param:
        arclength = point_1_param
    else:
        arclength = point_2_param
        
    intersection_x = x_dist + arclength * x_movement
    intersection_y = y_dist + arclength * y_movement
    
    # Recalculating, it's fine
    x_vel = (particle_2.velocity[0] - particle.velocity[0])
    y_vel = (particle_2.velocity[1] - particle.velocity[1])
    
    velocity_coefficient = 2 / ((particle.mass + particle_2.mass) * (intersection_x**2 + intersection_y**2)) * (x_vel * intersection_x + y_vel * intersection_y)
    v1_x = velocity_coefficient * particle_2.mass * intersection_x
    v1_y = velocity_coefficient * particle_2.mass * intersection_y
    v2_x = x_vel - velocity_coefficient * particle.mass * intersection_x
    v2_y = y_vel - velocity_coefficient * particle.mass * intersection_y
    
    # Move the particles to the collision point
    arclength_scaled = arclength * timestep_size
    num_grid_cols, num_grid_rows, arraygrid = particle.move_particle([particle.position[0] + arclength_scaled * particle.velocity[0], particle.position[1] + arclength_scaled * particle.velocity[1]],
        grid_cell_height, grid_cell_width, num_grid_cols, num_grid_rows, arraygrid)
    num_grid_cols, num_grid_rows, arraygrid = particle_2.move_particle([particle_2.position[0] + arclength_scaled * particle_2.velocity[0], particle_2.position[1] + arclength_scaled * particle_2.velocity[1]],
        grid_cell_height, grid_cell_width, num_grid_cols, num_grid_rows, arraygrid)
    
    # Update velocities, converting back to absolute velocities
    particle_2.velocity[0] = particle.velocity[0] + v2_x
    particle_2.velocity[1] = particle.velocity[1] + v2_y
    
    particle.velocity[0] = particle.velocity[0] + v1_x
    particle.velocity[1] = particle.velocity[1] + v1_y
    
    return arraygrid, Collided.COLLISION, particle_2, arclength, num_grid_cols, num_grid_rows


def get_arclength_to_wall(particle, wall, timestep_size):
    
    def get_intercept_coeff(particle, wall_line_vector, wall):
        
        # We want to solve a + tb = x + ky
        # a1 + tb1 - x1 = ky1
        # a2 + tb2 - x2 = ky2
        # t = 1/b1 * (x1 + ky1 - a1)
        # ky2 = a2 - x2 + b2/b1 * (x1 + ky1 - a1)
        # k(y2 - b2/b1y1) = a2 - x2 + b2/b1 * (x1 - a1)
        
        # Division by zero avoidance
        if particle.velocity[0] == 0 and particle.velocity[1] == 0:
            # Give up
            intercept_coeff = timestep_size + 1
        elif wall_line_vector[0] == 0:
            intercept_coeff = (particle.position[0] - wall[0]) / (-particle.velocity[0])
        else:
            wall_line_vector_ratio = wall_line_vector[1] / wall_line_vector[0]
            if particle.velocity[1] - wall_line_vector_ratio * particle.velocity[0] == 0:
                intercept_coeff = timestep_size + 1
            else:
                intercept_coeff = (wall[1] - particle.position[1] + wall_line_vector_ratio * (particle.position[0] - wall[0])) / ((particle.velocity[1] - wall_line_vector_ratio * particle.velocity[0]))
        return intercept_coeff
    
    wall_line_vector = [wall[2] - wall[0], wall[3] - wall[1]]
    magnitude = (wall_line_vector[0] ** 2 + wall_line_vector[1] ** 2) ** 0.5
    wall_line_vector = [i/magnitude for i in wall_line_vector]
    
    
    # I'm so glad I'm doing this in 2D
    normal_1 = [wall_line_vector[1], -wall_line_vector[0]]
    inflated_wall_1 = [wall[0] + particle.radius * normal_1[0], wall[1] + particle.radius * normal_1[1], wall[2] + particle.radius * normal_1[0], wall[3] + particle.radius * normal_1[1]]
    intercept_coeff_1 = get_intercept_coeff(particle, wall_line_vector, inflated_wall_1)
    
    normal_2 = [-wall_line_vector[1], wall_line_vector[0]]
    inflated_wall_2 = [wall[0] + particle.radius * normal_2[0], wall[1] + particle.radius * normal_2[1], wall[2] + particle.radius * normal_2[0], wall[3] + particle.radius * normal_2[1]]
    intercept_coeff_2 = get_intercept_coeff(particle, wall_line_vector, inflated_wall_2)
    
    # Technically this interacts incorrectly with the corners of walls, but like, I'm good. If the ball is going straight along a wall, it deserves it.
    
    return intercept_coeff_1, intercept_coeff_2, normal_1, normal_2

def check_for_wall_collisions_and_update(particle, walls, movement_remaining, num_grid_cols, num_grid_rows, arraygrid, timestep_size):
    # Note that this requires you to input all walls at once, since it mutates the particles
    
    intercept_coeffs = []
    intercepted_walls = []
    intercepted_wall_normals = []
    
    time_remaining = timestep_size * movement_remaining
    
    for wall in walls:
        wall_line_vector = [wall[2] - wall[0], wall[3] - wall[1]]
        magnitude = (wall_line_vector[0] ** 2 + wall_line_vector[1] ** 2) ** 0.5
        wall_line_vector = [i/magnitude for i in wall_line_vector]
        
        intercept_coeff_1, intercept_coeff_2, normal_1, normal_2 = get_arclength_to_wall(particle, wall, timestep_size)
        
        # Not a strict inequality, or else things will get stuck at walls
        # If we have to go backwards, or have to go forwards too much, we're done
        if intercept_coeff_1 <= 0 and intercept_coeff_2 <= 0 or intercept_coeff_1 > time_remaining and intercept_coeff_2 > time_remaining:
            continue
        
        elif intercept_coeff_2 <= 0:
            intercept_coeffs.append(intercept_coeff_1)
            intercepted_wall_normals.append(normal_1)
        
        elif intercept_coeff_1 < intercept_coeff_2:
            intercept_coeffs.append(intercept_coeff_1)
            intercepted_wall_normals.append(normal_1)
            
        else:
            intercept_coeffs.append(intercept_coeff_2)
            intercepted_wall_normals.append(normal_2)
        
        intercepted_walls.append(wall)
        
    
    nearest_intercept_coeff = time_remaining
    nearest_wall_normal = None
    nearest_wall = None
    
    for i in range(len(intercept_coeffs)):
        # Allowing for equality is important here
        if intercept_coeffs[i] <= nearest_intercept_coeff:
            nearest_intercept_coeff = intercept_coeffs[i]
            nearest_wall_normal = intercepted_wall_normals[i]
            nearest_wall = intercepted_walls[i]
    
    if nearest_wall is None:
        return Collided.NO_COLLISION, time_remaining
    
    
    hitpos = [particle.position[0] + particle.velocity[0] * nearest_intercept_coeff, particle.position[1] + particle.velocity[1] * nearest_intercept_coeff]
    
    # Scuffed check for if it's moveable or not, to account for relative velocities
    if len(nearest_wall) > 4:
        velocity_dot = (particle.velocity[0] - nearest_wall[5]) * nearest_wall_normal[0] + (particle.velocity[1] - nearest_wall[6]) * nearest_wall_normal[1]
    else:
        velocity_dot = particle.velocity[0] * nearest_wall_normal[0] + particle.velocity[1] * nearest_wall_normal[1]
    
    # Note that the normal by definition faces outwards, so the dot will be negative, and then flips again when multiplied by the nearest wall normal
    
    dvx = 2 * velocity_dot * nearest_wall_normal[0]
    dvy = 2 * velocity_dot * nearest_wall_normal[1]
    
    particle.velocity[0] = particle.velocity[0] - dvx
    particle.velocity[1] = particle.velocity[1] - dvy
    
    num_grid_cols, num_grid_rows, arraygrid = particle.move_particle(hitpos, grid_cell_height, grid_cell_width, num_grid_cols, num_grid_rows, arraygrid)
    
    # If the wall is stationary this impulse is just ignored
    impulse = [particle.mass * dvx, particle.mass * dvy]
    return Collided.COLLISION, 0, impulse, nearest_wall, num_grid_cols, num_grid_rows, arraygrid
    
    

def run_iteration(particles, num_grid_cols, num_grid_rows, arraygrid, timestep_size):
    
    particles_already_interacted = []
    arclengths_remaining = []
    
    for moveable_wall in moveable_walls:
        x_vel = moveable_wall[5]
        y_vel = moveable_wall[6]
        moveable_wall[0] = moveable_wall[0] + x_vel * timestep_size
        moveable_wall[1] = moveable_wall[1] + y_vel * timestep_size
        moveable_wall[2] = moveable_wall[2] + x_vel * timestep_size
        moveable_wall[3] = moveable_wall[3] + y_vel * timestep_size
    
    for particle in particles:
        
        # Since collisions will do an early evaluation of a particle's motion
        if particle in particles_already_interacted:
            continue
        
        
        cells_to_check = get_cells_to_check(particle.cell_coords[0], particle.cell_coords[1], arraygrid, *grid_bounds)
        
        # Hardcode checking the current cell first, quite ugly
        # Init array this way to avoid mutating the arraygrid. Could also copy via slicing.
        particles_to_check = []
        particles_to_check = particles_to_check + arraygrid[particle.cell_coords[0] - grid_bounds[0]][particle.cell_coords[1] - grid_bounds[1]]
        for cell_index in cells_to_check:
            # Don't consider the home cell twice
            if cell_index[0] == particle.cell_coords[0] - grid_bounds[0] and cell_index[1] == particle.cell_coords[1] - grid_bounds[1]:
                continue
            particles_to_check = particles_to_check + arraygrid[cell_index[0] - grid_bounds[0]][cell_index[1] - grid_bounds[1]]
        
        output = check_for_double_particle_collisions(particle, particles_to_check, num_grid_cols, num_grid_rows, arraygrid, timestep_size)
        
        arraygrid = output[0]
        particles_already_interacted.append(particle)
        if output[1] == Collided.COLLISION:
            particles_already_interacted.append(output[2])
            arclengths_remaining = arclengths_remaining + [1 - output[3], 1 - output[3]]
            
            num_grid_cols, num_grid_rows = output[4], output[5]
        
        else:
            arclengths_remaining.append(1)
        
        # After checking for collisions with particles, check for collisions with walls. Make the approximation that a particle
        # stops moving that frame upon collision with a wall. This will make my life a lot easier
    
    for i in range(len(particles_already_interacted)):
        output = check_for_wall_collisions_and_update(particles_already_interacted[i], walls + moveable_walls, arclengths_remaining[i], num_grid_cols, num_grid_rows, arraygrid, timestep_size)
        if output[0] == Collided.COLLISION:
            
            num_grid_cols, num_grid_rows, arraygrid = output[4], output[5], output[6]
            
            if output[3] in moveable_walls:
                # Update wall velocity by impulse over wall mass
                # I know this indexing is a bit horrifying but bear with me
                output[3][5] = output[3][5] + output[2][0] / output[3][4]
                output[3][6] = output[3][6] + output[2][1] / output[3][4]
        
        else:
            # ALl this work, just to hit this case most of the time, huh
            new_pos = [particles_already_interacted[i].position[0] + particles_already_interacted[i].velocity[0] * output[1],
                particles_already_interacted[i].position[1] + particles_already_interacted[i].velocity[1] * output[1]]
            
            num_grid_cols, num_grid_rows, arraygrid = particles_already_interacted[i].move_particle(new_pos, grid_cell_height, grid_cell_width, num_grid_cols, num_grid_rows, arraygrid)
    
    return num_grid_cols, num_grid_rows, arraygrid


particle_1 = Particle([1, 0], [1000000, 0], 1, 1, grid_cell_height, grid_cell_width)
particle_2 = Particle([3, 10], [0, 0], 1, 1, grid_cell_height, grid_cell_width)

arraygrid[particle_1.cell_coords[0] - grid_bounds[0]][particle_1.cell_coords[1] - grid_bounds[1]].append(particle_1)
arraygrid[particle_2.cell_coords[0] - grid_bounds[0]][particle_2.cell_coords[1] - grid_bounds[1]].append(particle_2)

particles = [particle_1, particle_2]

for i in range(100):
    num_grid_cols, num_grid_rows, arraygrid = run_iteration(particles, num_grid_cols, num_grid_rows, arraygrid, 0.01)
    print(len(arraygrid))
    print("Particle 1 position is", particle_1.position)
    print("Particle 1 velocity is", particle_1.velocity)
    print(moveable_walls)

# TODO swap the moveable wall to the left and make it work

#print(check_wall(Particle([0, 5], [5, -2], 1, 1, grid_cell_height, grid_cell_width), [0, 0, 10, 0], 1))