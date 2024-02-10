# Define the parameters for the sphere
radius <- 1  # Radius of the sphere
resolution <- 50  # Number of vertices to use for the sphere

# Generate the coordinates for the sphere
theta <- seq(0, 2*pi, length.out = resolution)
phi <- seq(0, pi, length.out = resolution)
coords <- expand.grid(theta = theta, phi = phi)

# Convert spherical coordinates to Cartesian coordinates
x <- radius * sin(coords$phi) * cos(coords$theta)
y <- radius * sin(coords$phi) * sin(coords$theta)
z <- radius * cos(coords$phi)

sphere_df <- tibble::tibble(x1 = x, x2 = y, x3 = z)

langevitour::langevitour(sphere_df)

# library(rgl)
# open3d()
# shade3d(
#   spheres3d(
#     x = x,
#     y = y,
#     z = z,
#     radius = 0.05,  # Adjust the size of the spheres for visualization
#     col = "blue"   # Adjust the color of the spheres
#   )
# )

# # Install and load the pracma package
# install.packages("pracma")
# library(pracma)
#
# # Define the radius of the sphere
# radius <- 1
#
# # Generate coordinates for the sphere
# theta <- seq(0, pi, length.out = 30)
# phi <- seq(0, 2 * pi, length.out = 30)
#
# # Generate grid of theta and phi values
# grid <- expand.grid(theta = theta, phi = phi)
#
# # Convert spherical coordinates to Cartesian coordinates
# sphere_coords <- cbind(
#   x = radius * sin(grid$theta) * cos(grid$phi),
#   y = radius * sin(grid$theta) * sin(grid$phi),
#   z = radius * cos(grid$theta)
# )
#
# langevitour::langevitour(sphere_coords)
