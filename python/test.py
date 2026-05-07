import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("assets/world_traj.csv")

earth = df[["earth_x", "earth_y", "earth_z"]].values
urath = df[["urath_x", "urath_y", "urath_z"]].values
sat = df[["sat_x", "sat_y", "sat_z"]].values

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection="3d")

ax.plot(earth[:, 0], earth[:, 1], earth[:, 2], label="Earth", color="blue")
ax.plot(urath[:, 0], urath[:, 1], urath[:, 2], label="Urath", color="red")
ax.plot(sat[:, 0], sat[:, 1], sat[:, 2], label="Satellite", color="green")

all_points = df[
    [
        "earth_x",
        "earth_y",
        "earth_z",
        "urath_x",
        "urath_y",
        "urath_z",
        "sat_x",
        "sat_y",
        "sat_z",
    ]
].values

max_range = (all_points.max(axis=0) - all_points.min(axis=0)).max() / 2.0
mid_x = (all_points[:, 0].max() + all_points[:, 0].min()) * 0.5
mid_y = (all_points[:, 1].max() + all_points[:, 1].min()) * 0.5
mid_z = (all_points[:, 2].max() + all_points[:, 2].min()) * 0.5

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

ax.legend()
plt.tight_layout()
plt.show()
