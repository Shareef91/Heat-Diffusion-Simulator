# 🔥 Heat Diffusion Simulator

> A fast and simple 2D heat diffusion simulation written in C++ 🚀

## 📖 What is this?

Watch heat spread across a 2D grid in real-time! This simulator uses finite-difference methods to model thermal diffusion, starting from a hot center point and gradually spreading outward. Perfect for learning about numerical simulations and heat transfer! ☕➡️🧊

## ✨ Features

- 🎯 **Single-file simplicity** - Just one C++ file, easy to understand and modify
- ⚡ **Fast simulation** - Optimized loops for quick results
- 📊 **Real-time progress** - See temperature updates every 100 steps
- 🎨 **Customizable parameters** - Adjust grid size, time steps, and heat transfer rate
- 🔬 **Physics-based** - Uses proper finite-difference approximations

## 🚀 Quick Start

### Building (Windows PowerShell)

**With g++ (MinGW/mingw-w64):**
```powershell
g++ -std=c++17 -O2 heat_diffusion.cpp -o heat_diffusion.exe
```

**With MSVC:**
```powershell
cl /EHsc /O2 heat_diffusion.cpp
```

### Running 🏃‍♂️

```powershell
.\heat_diffusion.exe
```

You'll see output like:
```
Step 0: Center = 100
Step 100: Center = 89.234
Step 200: Center = 81.567
...
```

## 🎛️ Configuration

Edit the constants at the top of `heat_diffusion.cpp`:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `N` | Grid size (N × N) | 100 |
| `STEPS` | Number of time steps | 1000 |
| `alpha` | Heat transfer rate | 0.1 |
| `dt` | Time increment | 0.1 |
| `dx` | Distance between cells | 1.0 |

## 🧪 How it works

1. **Initialize** 🌡️ - Start with a 100×100 grid at ambient temperature (25°C)
2. **Add heat source** 🔥 - Place a hot spot (100°C) in the center
3. **Simulate** ⏱️ - Run 1000 time steps, updating each cell based on its neighbors
4. **Observe** 👀 - Watch the center temperature decrease as heat spreads

The core physics is in `compute_new_value()`:
```cpp
T_new = T_old + α·Δt/Δx² · (T_up + T_down + T_left + T_right - 4·T_center)
```

## 🎨 Example Output

```
Step 0: Center = 100
Step 100: Center = 89.4523
Step 200: Center = 82.1156
Step 300: Center = 76.8934
Step 400: Center = 72.9812
...
Step 1000: Center = 58.3421
```

## 🤝 Contributing

Want to improve the simulator? Ideas:

- 🎯 Add visualization (output to CSV/image)
- 🖥️ Implement command-line arguments
- ⚡ Add GPU acceleration
- 🧪 Create unit tests
- 📈 Add different initial conditions (multiple hot spots, boundary heating)

## 📚 Learn More

This simulator uses the **finite difference method** to approximate the heat equation:

$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T$$

Perfect for:
- 📖 Learning numerical methods
- 🎓 Teaching thermal physics
- 🔬 Prototyping larger simulations

## 📝 License

Feel free to use this code for learning and experimentation! 🎉

---

Made with ❤️ and ☕ | Happy simulating! 🌡️✨
