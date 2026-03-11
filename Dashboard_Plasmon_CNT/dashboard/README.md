# Surface Polariton Dispersion Simulator (Web Dashboard)

An interactive, high-performance web dashboard for visual exploration of surface polariton dispersion and field profiles.

## 🚀 Getting Started

Launch the dashboard instantly by opening `index.html` in any modern web browser (Chrome, Firefox, Safari, or Edge).

```bash
# Example: Open using Chrome via terminal
start chrome index.html
```

## 🛠️ Key Features

### Interactive Controls
- **Material Parameters**: Adjust $\varepsilon_1$, $\varepsilon_3$, $\omega_p$, and $\Gamma$ using real-time sliders.
- **Film Geometry**: Dynamically modify the metal film thickness ($h$).
- **Frequency Picker**: Select a specific frequency to obtain exact numerical roots and physical metrics.

### Live Visualizations
- **Dispersion Plot**: Real-time graph of the propagation constant $\beta_R$ vs frequency.
- **Attenuation Plot**: Visualization of $\beta_I$ (losses) across the spectrum.
- **Field Profile**: Depth-resolved intensity profile $|f(z)|^2$ at the selected frequency.
- **Analytic Overlays**: Reference lines for single-interface limits and light lines.

## 🔬 Core Technologies
- **Vanilla JavaScript**: High-efficiency implementation without external dependencies.
- **Canvas API**: Custom plotting engine for smooth, reactive rendering.
- **Muller's Method**: Complex root-finding algorithm mirroring the C++ core.
- **Adaptive Layout**: Responsive design with support for Light and Dark modes.

## 📖 In-App Help
The dashboard includes built-in modals (Guide, Parameters, Curves, Algorithms, and Wiki) to provide on-the-fly scientific context and physical guidance.

---
*Created as part of the Surface Polariton Implementation project at UEF.*
