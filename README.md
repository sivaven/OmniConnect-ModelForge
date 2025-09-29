# OmniConnect-ModelForge

This repository is an **umbrella project** for computational neuroscience tools and models.  
It is designed to integrate multiple approaches to studying brain networks across scales, from **directed connectomes** to **spiking neuron models**.  

---

## Structure

All projects live under the `projects/` directory. Each project has its own README and its own `requirements.txt`.

```
OmniConnect-ModelForge/
├── README.md              # Top-level overview (this file)
├── LICENSE                # MIT license for the whole repository
├── requirements.txt       # Shared dependencies 
│
└── projects/
    ├── cross_species_connectomics/   # Directed connectomics across species
    └── project2/               
```

---

## Projects

### 1. Cross-Species Connectomics
**Path:** `projects/cross_species_connectomics/`  

Directed connectome framework integrating **mouse viral tracing** with **diffusion MRI** across mouse, marmoset, rhesus macaque, and human.  
Includes tractography shell scripts, connectome construction, statistical analyses, and figure generation from publication.  

📄 [Project README](projects/cross_species_connectomics/README.md)  
📖 Citation:  
> **Cross-species brain circuitry from diffusion MRI tractography and mouse viral tracing.**  
> Siva Venkadesh, Wen-Jieh Linn, Yuhe Tian, G Allan Johnson, Fang-Cheng Yeh.  
> *bioRxiv* (2025). doi: [10.1101/2025.09.07.674762](https://doi.org/10.1101/2025.09.07.674762)  

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<your-username>/OmniConnect-ModelForge.git
cd OmniConnect-ModelForge
```

Install shared dependencies:

```bash
pip install -r requirements.txt
```

Each project may also define its own `requirements.txt` for additional dependencies.

---

## License

This repository is released under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.

---
