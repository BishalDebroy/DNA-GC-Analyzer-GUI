#!/usr/bin/env python3
"""
DNA Sequence Analyzer – GUI with embedded plot (Tkinter + Matplotlib)
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import Counter
import matplotlib
matplotlib.use("TkAgg")                     # Force Tkinter backend (no Qt needed)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ----------------------------------------------------------------------
# Core analysis (unchanged from CLI version)
# ----------------------------------------------------------------------
def parse_fasta(text):
    """Parse FASTA text (str) and return list of (header, seq) tuples."""
    records = []
    header = None
    seq_parts = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts).upper()))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records

def compute_stats(header, sequence):
    total_len = len(sequence)
    counts = Counter(sequence)
    a = counts.get("A", 0)
    t = counts.get("T", 0)
    g = counts.get("G", 0)
    c = counts.get("C", 0)
    other = total_len - (a + t + g + c)
    valid = a + t + g + c
    gc = (g + c) / valid * 100 if valid else 0.0
    at = (a + t) / valid * 100 if valid else 0.0
    return {
        "header": header,
        "length": total_len,
        "A": a, "T": t, "G": g, "C": c,
        "other": other,
        "GC": gc,
        "AT": at
    }

# ----------------------------------------------------------------------
# GUI application
# ----------------------------------------------------------------------
class DNAAnalyzerApp:
    def __init__(self, root):
        self.root = root
        root.title("DNA GC Analyzer & Comparator")
        root.geometry("950x700")
        root.minsize(800, 600)

        # Variables
        self.stats_list = []          # holds dicts after analysis
        self.fig = None
        self.canvas = None

        # --- Menu bar ---
        menubar = tk.Menu(root)
        root.config(menu=menubar)
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open FASTA...", command=self.load_file)
        file_menu.add_command(label="Export Plot...", command=self.export_plot)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=root.quit)
        menubar.add_cascade(label="File", menu=file_menu)

        # --- Input area ---
        input_frame = ttk.LabelFrame(root, text="Input (paste or load a multi-FASTA)")
        input_frame.pack(fill=tk.BOTH, expand=False, padx=10, pady=(10, 5))

        self.text_input = tk.Text(input_frame, height=8, wrap=tk.NONE)
        self.text_input.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        scroll_y = ttk.Scrollbar(input_frame, orient=tk.VERTICAL, command=self.text_input.yview)
        scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.text_input.configure(yscrollcommand=scroll_y.set)

        # --- Buttons ---
        btn_frame = ttk.Frame(root)
        btn_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Button(btn_frame, text="Analyze", command=self.run_analysis).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear", command=self.clear_all).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Load Example", command=self.load_example).pack(side=tk.LEFT, padx=5)

        # --- Results table ---
        table_frame = ttk.LabelFrame(root, text="Sequence statistics")
        table_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 5))

        columns = ("Entry ID", "Length", "A", "T", "G", "C", "Other", "GC%", "AT%")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=6)
        for col in columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=90, anchor=tk.CENTER)
        self.tree.column("Entry ID", width=200, anchor=tk.W)   # wider for names

        vsb = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        table_frame.grid_rowconfigure(0, weight=1)
        table_frame.grid_columnconfigure(0, weight=1)

        # --- Plot area ---
        plot_frame = ttk.LabelFrame(root, text="GC vs AT comparison")
        plot_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))

        self.plot_container = ttk.Frame(plot_frame)
        self.plot_container.pack(fill=tk.BOTH, expand=True)

        # Status bar
        self.status = ttk.Label(root, text="Ready. Paste or load a FASTA file and click Analyze.",
                                relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(fill=tk.X, side=tk.BOTTOM, ipady=2)

    # ------------------------------------------------------------------
    def load_file(self):
        """Open a FASTA file and insert its content into the text widget."""
        filepath = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")]
        )
        if not filepath:
            return
        try:
            with open(filepath, "r") as f:
                content = f.read()
            self.text_input.delete("1.0", tk.END)
            self.text_input.insert("1.0", content)
            self.status.config(text=f"Loaded {filepath}")
        except Exception as e:
            messagebox.showerror("Error", f"Could not read file:\n{e}")

    def load_example(self):
        """Insert a short demo multi-FASTA."""
        example = """>Escherichia coli
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
>Homo sapiens
CCTAGCTAGCTAGCTAATATATATGCGCGCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGC
>Saccharomyces cerevisiae
AAAAAAAAATTTTTTTTTTGGGGGGGGGCCCCCCCCC
"""
        self.text_input.delete("1.0", tk.END)
        self.text_input.insert("1.0", example)
        self.status.config(text="Example loaded. Click Analyze.")

    def clear_all(self):
        """Clear input, table, and plot."""
        self.text_input.delete("1.0", tk.END)
        self.tree.delete(*self.tree.get_children())
        self.clear_plot()
        self.stats_list = []
        self.status.config(text="Cleared.")

    def clear_plot(self):
        """Remove existing plot from the GUI."""
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        plt.close("all")

    def run_analysis(self):
        """Parse input, compute stats, populate table, and draw plot."""
        fasta_text = self.text_input.get("1.0", tk.END).strip()
        if not fasta_text:
            messagebox.showwarning("No input", "Please paste or load a FASTA sequence.")
            return

        records = parse_fasta(fasta_text)
        if not records:
            messagebox.showwarning("No sequences", "No valid FASTA records found.")
            return

        self.stats_list = [compute_stats(h, s) for h, s in records]

        # Update table
        self.tree.delete(*self.tree.get_children())
        for s in self.stats_list:
            self.tree.insert("", tk.END, values=(
                s["header"], s["length"],
                s["A"], s["T"], s["G"], s["C"],
                s["other"],
                f"{s['GC']:.2f}", f"{s['AT']:.2f}"
            ))

        # Redraw plot
        self.draw_plot()
        self.status.config(text=f"Analyzed {len(self.stats_list)} sequence(s).")

    def draw_plot(self):
        """Create (or recreate) the Matplotlib bar chart inside the plot frame."""
        self.clear_plot()   # remove old canvas

        if not self.stats_list:
            return

        # Prepare data
        labels = [s["header"].split()[0] for s in self.stats_list]   # first word
        gc_vals = [s["GC"] for s in self.stats_list]
        at_vals = [s["AT"] for s in self.stats_list]

        # Create figure and axes
        self.fig, ax = plt.subplots(figsize=(6, 3.5))   # size in inches
        x = range(len(labels))
        width = 0.35

        bars_gc = ax.bar([i - width/2 for i in x], gc_vals, width, label="GC%", color="#2196F3")
        bars_at = ax.bar([i + width/2 for i in x], at_vals, width, label="AT%", color="#FF9800")

        # Annotate bars
        for bar in bars_gc:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 0.5, f"{h:.1f}",
                    ha="center", va="bottom", fontsize=8)
        for bar in bars_at:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., h + 0.5, f"{h:.1f}",
                    ha="center", va="bottom", fontsize=8)

        ax.set_ylabel("Percentage (%)")
        ax.set_title("GC and AT Content Comparison")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.legend()
        ax.set_ylim(0, 105)
        ax.grid(axis="y", linestyle="--", alpha=0.7)
        plt.tight_layout()

        # Embed in Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_container)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def export_plot(self):
        """Save the current plot as a PNG file."""
        if not self.fig:
            messagebox.showinfo("No plot", "Analyze some sequences first to generate a plot.")
            return
        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("All files", "*.*")]
        )
        if filepath:
            self.fig.savefig(filepath, dpi=150)
            self.status.config(text=f"Plot saved to {filepath}")
            messagebox.showinfo("Saved", f"Plot exported as {filepath}")

# ----------------------------------------------------------------------
# Main entry point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    app = DNAAnalyzerApp(root)
    root.mainloop()