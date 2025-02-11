import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime
import matplotlib as mpl
mpl.use('Qt5Agg')
# Simulation parameters
start_date = datetime.datetime(2013, 11, 4)
time_step_days = 1

# Sample data (Replace with real data)
n = 1000  # Number of time steps
t = np.arange(n)
positions_asteroid = np.random.rand(3, n) * 5  # Fake asteroid data
positions_planets = np.random.rand(12, n) * 5  # Fake planets data

# Extract individual planet positions
earth_positions = positions_planets[0:3, :]
mars_positions = positions_planets[3:6, :]
jupiter_positions = positions_planets[6:9, :]
saturn_positions = positions_planets[9:12, :]

# Initialize figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_xlabel("X (AU)")
ax.set_ylabel("Y (AU)")

# Colors
colors = {
    "earth": "#3b7ddd",  # Blue
    "mars": "#c1440e",   # Red-orange
    "jupiter": "#daa520",  # Golden-brown
    "saturn": "#f4a460",  # Sandy yellow
    "asteroid": "#555555"  # Dark gray
}

# Plot elements
asteroid_trail, = ax.plot([], [], '-', color=colors["asteroid"], alpha=0.6, label="Asteroid 2007 W266")
earth_trail, = ax.plot([], [], '-', color=colors["earth"], alpha=0.6, label="Earth")
mars_trail, = ax.plot([], [], '-', color=colors["mars"], alpha=0.6, label="Mars")
jupiter_trail, = ax.plot([], [], '-', color=colors["jupiter"], alpha=0.6, label="Jupiter")
saturn_trail, = ax.plot([], [], '-', color=colors["saturn"], alpha=0.6, label="Saturn")

asteroid_dot, = ax.plot([], [], 'o', color=colors["asteroid"], markersize=6)
earth_dot, = ax.plot([], [], 'o', color=colors["earth"], markersize=6)
mars_dot, = ax.plot([], [], 'o', color=colors["mars"], markersize=6)
jupiter_dot, = ax.plot([], [], 'o', color=colors["jupiter"], markersize=6)
saturn_dot, = ax.plot([], [], 'o', color=colors["saturn"], markersize=6)
sun_dot = ax.scatter(0,0,color="y",marker="o", s=60, label="Sun")
# Title
title = ax.set_title("")

# Update function for animation
def update(frame):
    # Update trails
    asteroid_trail.set_data(positions_asteroid[0, :frame], positions_asteroid[1, :frame])
    earth_trail.set_data(earth_positions[0, :frame], earth_positions[1, :frame])
    mars_trail.set_data(mars_positions[0, :frame], mars_positions[1, :frame])
    jupiter_trail.set_data(jupiter_positions[0, :frame], jupiter_positions[1, :frame])
    saturn_trail.set_data(saturn_positions[0, :frame], saturn_positions[1, :frame])

    # Update moving points
    asteroid_dot.set_data(positions_asteroid[0, frame], positions_asteroid[1, frame])
    earth_dot.set_data(earth_positions[0, frame], earth_positions[1, frame])
    mars_dot.set_data(mars_positions[0, frame], mars_positions[1, frame])
    jupiter_dot.set_data(jupiter_positions[0, frame], jupiter_positions[1, frame])
    saturn_dot.set_data(saturn_positions[0, frame], saturn_positions[1, frame])

    # Update title with month and year
    current_date = start_date + datetime.timedelta(days=frame * time_step_days)
    title.set_text(current_date.strftime("%B %Y"))

    return (asteroid_trail, earth_trail, mars_trail, jupiter_trail, saturn_trail,
            asteroid_dot, earth_dot, mars_dot, jupiter_dot, saturn_dot, title)

# Create animation
ani = animation.FuncAnimation(fig, update, frames=n, interval=2, repeat=True)

# Show animation
plt.legend()
plt.grid()
plt.show()
