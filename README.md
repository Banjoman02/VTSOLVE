# VT-SOLVE

**Virginia Tech Space Object Location and Velocity Estimator**

VT-SOLVE is a Python-based orbit determination tool developed at Virginia Tech. It is designed to estimate the orbits of asteroids and other space objects orbiting the Sun using optical images captured by ground-based telescopes.

---

## 🚀 Overview

VT-SOLVE performs orbit determination by analyzing a series of images containing the target object. The software implements a platesolving pipeline to determine accurate celestial coordinates (RA/Dec) of the object in each image. Using three such observations, it then applies **Gaussian Initial Orbit Determination (IOD)** to estimate the object's orbital elements around the Sun.

---

## 🔍 Key Features

- 📷 **Platesolving**: Convert pixel coordinates of the object in astronomical images into celestial coordinates using advanced platesolving techniques.
- 🌌 **Asteroid Tracking**: Detect and track asteroids across multiple images.
- 🪐 **Gaussian IOD**: Compute orbital elements from three distinct observations using the classical Gaussian method.
- 🛰️ **Solar Orbit Focus**: Tailored specifically for orbit determination of objects around the Sun.

---

## System Requirements.

VT-SOLVE is supported on Linux systems with Python 3.9 or later installed. It is not currently supported on Windows or MacOS.

## 🛠 Installation

1. **Create a Python Virtual Environment** (recommended):

```bash
python -m venv venv
source venv/bin/activate
```

2. **Install VTSOLVE**

```bash
make install
```

## Usage

No entry point and usage exists at the moment.