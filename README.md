# DNA GC Content Analyzer (GUI)

A Python-based bioinformatics tool that accepts DNA sequences in **FASTA format** and computes:
- Sequence length
- Nucleotide frequencies (A, T, G, C, other)
- **GC content** (%)
- **AT content** (%)

Results are displayed in a **sortable table** and a **bar chart** comparing GC% and AT% across all entries.  
Built with **Tkinter** and **Matplotlib** – no Qt dependencies required.

---

## Features

- Paste multi‑FASTA directly or load from a file
- Instant statistics for every entry
- Embedded **grouped bar chart** with value labels
- Export plots as **PNG** or **PDF**
- Clean, cross‑platform graphical interface

---

## Screenshots

![GUI screenshot](screenshot.png)   <!-- Replace with actual screenshot if available -->

---

## Requirements

- Python 3.6 or higher
- **Matplotlib** (install via pip)
- **Tkinter** – included with most Python distributions.  
  if not, on **Ubuntu/Debian**, install it with:
  ```bash
  sudo apt update && sudo apt install python3-tk
  ```
## Installation
1. Clone the repository in your linux
   ```
   git clone https://github.com/BishalDebroy/DNA-GC-Analyzer-GUI.git
   cd DNA-GC-Analyzer-GUI
   ```

 2. Create a virtual enveronment (recomanded)
    ```
    conda create -n venv
    ```
  3. Install matplotlib v3.5.0
     ```
     pip install matplotlib>==3.5.0
     ```
   4. Run the program
      ```
      python3 dna_gc_analyzer_gui.py
      ```
  5. to kill the program, in the terminal press **ctrl+c**
