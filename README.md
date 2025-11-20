# Asteroid Rendezvous Launch Window Tool

This project computes the **optimal launch date** from Earth to a target asteroid using classical orbital elements and a Lambert transfer.  
It is fully interactive — the user enters asteroid elements, mission dates, and flight time, and the script computes the Δv curve and best solution.

This version is compatible with **poliastro 0.7.0**, which is used automatically on macOS ARM + Python 3.13.

---

##  Features

- Build asteroid orbit from user-entered classical elements  
- Earth ephemeris from poliastro  
- Propagate asteroid to arrival date  
- Solve Lambert transfers for multiple departure dates  
- Compute:
  - Δv at departure  
  - Δv at arrival  
  - total Δv  
- Scan a launch window  
- Plot Δv vs departure date  
- Print the optimal launch date and arrival date  

---

## Example Inputs (used during development)

### **Asteroid orbital elements**
a = 1.8 AU
e = 0.2
i = 10 deg
RAAN Ω = 80 deg
ω = 30 deg
M = 150 deg
Epoch = 2025-01-01 00:00

### **Mission parameters**

Launch window start: 2026-01-01 00:00
Launch window end: 2026-12-31 00:00
Time of flight: 200 days
Number of samples: 30 departure dates

---

## Example Output


Best departure: 2026-08-27 11:35 TDB
Arrival: 2027-03-15 11:35 TDB
Minimum Δv: 37.174 km/s


A Δv vs departure-date plot is also generated.

---

## How to Run

### 1. Create a virtual environment
python3 -m venv .venv
source .venv/bin/activate

### 2. Install required packages
python -m pip install numpy matplotlib astropy poliastro

(Your system will automatically install **poliastro 0.7.0**, which matches the script.)

### 3. Run the program
(python asteroid_launch_window.py)



---

## 📁 Project Structure
asteroid_window/
│
├── asteroid_launch_window.py # main tool
└── README.md # this file


