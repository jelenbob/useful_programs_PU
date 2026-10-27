from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp


OUTPUT_DIR = Path(__file__).parent / "files"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

G = 9.81
V0 = 36.0
ALPHA_DEG = 7.0
MASS = 0.145
AIR_DENSITY = 1.225
DRAG_COEFFICIENT = 0.47
BALL_RADIUS = 0.0366
AREA = np.pi * BALL_RADIUS**2
DRAG_FACTOR = 0.5 * AIR_DENSITY * DRAG_COEFFICIENT * AREA

BLUE = "#2864a5"
ORANGE = "#d47722"
GREEN = "#168578"
PURPLE = "#8455a4"
INK = "#233044"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.2,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def simulate_with_drag(speed, angle_rad, sample_count=1000):
    initial = [0.0, 0.0, speed * np.cos(angle_rad), speed * np.sin(angle_rad)]

    def equations(_time, state):
        _x, _y, vx, vy = state
        velocity = np.hypot(vx, vy)
        return [
            vx,
            vy,
            -(DRAG_FACTOR / MASS) * velocity * vx,
            -G - (DRAG_FACTOR / MASS) * velocity * vy,
        ]

    def touchdown(_time, state):
        return state[1]

    touchdown.terminal = True
    touchdown.direction = -1
    trial = solve_ivp(
        equations,
        (0.0, 30.0),
        initial,
        events=touchdown,
        dense_output=True,
        max_step=0.025,
        rtol=1e-9,
        atol=1e-11,
    )
    if not trial.success or len(trial.t_events[0]) != 1:
        raise RuntimeError("Nie udało się wyznaczyć chwili upadku.")
    times = np.linspace(0.0, trial.t_events[0][0], sample_count)
    return times, trial.sol(times)


def save_figure(figure, name):
    figure.savefig(OUTPUT_DIR / name, dpi=180, bbox_inches="tight")
    figure.savefig(OUTPUT_DIR / name.replace(".png", ".svg"), bbox_inches="tight")
    plt.close(figure)


angle = np.deg2rad(ALPHA_DEG)
flight_time_ideal = 2 * V0 * np.sin(angle) / G
time_ideal = np.linspace(0.0, flight_time_ideal, 1000)
x_ideal = V0 * np.cos(angle) * time_ideal
y_ideal = V0 * np.sin(angle) * time_ideal - 0.5 * G * time_ideal**2
vx_ideal = np.full_like(time_ideal, V0 * np.cos(angle))
vy_ideal = V0 * np.sin(angle) - G * time_ideal

time_drag, state_drag = simulate_with_drag(V0, angle)
x_drag, y_drag, vx_drag, vy_drag = state_drag
speed_drag = np.hypot(vx_drag, vy_drag)
range_ideal = x_ideal[-1]
range_drag = x_drag[-1]
max_height_ideal = np.max(y_ideal)
max_height_drag = np.max(y_drag)

figure, axis = plt.subplots(figsize=(9, 5), layout="constrained")
axis.plot(x_ideal, y_ideal, color=BLUE, label="Bez oporu powietrza")
axis.plot(x_drag, y_drag, color=ORANGE, label="Z oporem kwadratowym")
axis.scatter([0, range_ideal, range_drag], [0, 0, 0], color=[INK, BLUE, ORANGE], s=24, zorder=3)
axis.set(xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]")
axis.set_title("Tor lotu: model idealny i opór powietrza", loc="left")
axis.legend()
save_figure(figure, "tor_lotu.png")

figure, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
axes[0].plot(time_ideal, x_ideal, color=BLUE, label="Bez oporu")
axes[0].plot(time_drag, x_drag, color=ORANGE, label="Z oporem")
axes[0].set_ylabel("x [m]")
axes[0].set_title("Współrzędne położenia w czasie", loc="left")
axes[0].legend()
axes[1].plot(time_ideal, y_ideal, color=BLUE, label="Bez oporu")
axes[1].plot(time_drag, y_drag, color=ORANGE, label="Z oporem")
axes[1].set(xlabel="Czas t [s]", ylabel="y [m]")
save_figure(figure, "polozenie_czas.png")

figure, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
axes[0].plot(time_ideal, vx_ideal, color=BLUE, label="vx bez oporu")
axes[0].plot(time_drag, vx_drag, color=ORANGE, label="vx z oporem")
axes[0].set_ylabel("Prędkość pozioma [m/s]")
axes[0].set_title("Składowe prędkości i szybkość", loc="left")
axes[0].legend(ncol=2)
axes[1].plot(time_ideal, vy_ideal, color=GREEN, label="vy bez oporu")
axes[1].plot(time_drag, vy_drag, color=PURPLE, label="vy z oporem")
axes[1].plot(time_drag, speed_drag, color=INK, linestyle="--", label="Szybkość z oporem")
axes[1].axhline(0, color=INK, linewidth=0.8)
axes[1].set(xlabel="Czas t [s]", ylabel="Prędkość [m/s]")
axes[1].legend(ncol=2)
save_figure(figure, "predkosc_czas.png")

kinetic_energy = 0.5 * MASS * speed_drag**2
potential_energy = MASS * G * y_drag
mechanical_energy = kinetic_energy + potential_energy
dissipated_energy = cumulative_trapezoid(
    DRAG_FACTOR * speed_drag**3, time_drag, initial=0.0
)
energy_residual = mechanical_energy + dissipated_energy - mechanical_energy[0]

figure, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout="constrained")
axes[0].plot(time_drag, kinetic_energy, color=BLUE, label="Energia kinetyczna")
axes[0].plot(time_drag, potential_energy, color=GREEN, label="Energia potencjalna grawitacji")
axes[0].plot(time_drag, mechanical_energy, color=ORANGE, label="Energia mechaniczna")
axes[0].set_ylabel("Energia [J]")
axes[0].set_title("Energia przy działaniu oporu powietrza", loc="left")
axes[0].legend(ncol=2)
axes[1].plot(time_drag, dissipated_energy, color=PURPLE, label="Praca utracona na opór")
axes[1].plot(time_drag, energy_residual, color=INK, linestyle="--", label="Reszta bilansu energii")
axes[1].set(xlabel="Czas t [s]", ylabel="Energia [J]")
axes[1].legend()
save_figure(figure, "bilans_energii.png")

angle_values = np.linspace(np.deg2rad(2.0), np.deg2rad(80.0), 157)
ranges_ideal = V0**2 * np.sin(2 * angle_values) / G
ranges_drag = np.array([
    simulate_with_drag(V0, value, sample_count=2)[1][0, -1]
    for value in angle_values
])
best_angle = angle_values[np.argmax(ranges_drag)]

figure, axis = plt.subplots(figsize=(9, 5), layout="constrained")
axis.plot(angle_values, ranges_ideal, color=BLUE, label="Bez oporu powietrza")
axis.plot(angle_values, ranges_drag, color=ORANGE, label="Z oporem kwadratowym")
axis.axvline(np.pi / 4, color=BLUE, linestyle=":", alpha=0.7,
             label="0,785 rad: maksimum idealne")
axis.axvline(best_angle, color=ORANGE, linestyle="--", alpha=0.7,
             label=f"{best_angle:.3f} rad: maksimum z oporem")
axis.scatter([angle], [range_drag], color=INK, zorder=3,
             label=f"Badany kąt {angle:.3f} rad")
axis.set(xlabel="Kąt wyrzutu α [rad]", ylabel="Zasięg [m]", xlim=(0, np.pi / 2))
axis.set_title("Zasięg w funkcji kąta wyrzutu", loc="left")
axis.legend(fontsize=9, ncol=2)
save_figure(figure, "zasieg_kat.png")

assert np.max(np.abs(energy_residual)) < 2e-3
assert range_drag < range_ideal
print(f"v0 = {V0:.2f} m/s; alpha = {angle:.4f} rad ({ALPHA_DEG:.2f} deg)")
print(f"Pole przekroju kuli = {AREA:.6f} m^2; czynnik oporu = {DRAG_FACTOR:.6f} kg/m")
print(f"Bez oporu: czas lotu = {flight_time_ideal:.4f} s, zasięg = {range_ideal:.4f} m, ymax = {max_height_ideal:.4f} m")
print(f"Z oporem: czas lotu = {time_drag[-1]:.4f} s, zasięg = {range_drag:.4f} m, ymax = {max_height_drag:.4f} m")
print(f"Kąt maksymalnego zasięgu z oporem = {best_angle:.4f} rad ({np.rad2deg(best_angle):.2f} deg)")
print(f"Maksymalna reszta bilansu energii = {np.max(np.abs(energy_residual)):.3e} J")
print(f"Zapisano wykresy w: {OUTPUT_DIR}")