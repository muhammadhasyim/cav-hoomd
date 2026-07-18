#!/usr/bin/env python3
"""
LaTeX font configuration for Adobe Illustrator compatibility.

This module provides a standardized LaTeX configuration that ensures 
proper Computer Modern font embedding in PDF files, making them compatible 
with Adobe Illustrator.

Key fixes:
1. Removes 'lmodern' package that creates LMMathItalic/LMMathSymbols 
2. Uses pure Computer Modern fonts
3. Configures Type 1 font embedding (not Type 3)
4. Proper amsmath packages for mathematical symbols

Usage:
    from latex_config_adobe import setup_latex_fonts, latex_safe
    
    use_latex = setup_latex_fonts()
    plt.title(latex_safe(r'$\\tau = 2\\pi f$', 'τ = 2πf', use_latex))
"""

import os


def _sanitize_tex_environment():
    """Conda's LD_LIBRARY_PATH breaks luatex/kpathsea used by matplotlib usetex."""
    os.environ.pop("LD_LIBRARY_PATH", None)
    os.environ.pop("LD_PRELOAD", None)


def _reset_matplotlib_tex_caches():
    """Clear cached (possibly broken) luatex kpathsea subprocesses."""
    try:
        from matplotlib import dviread
        dviread._LuatexKpsewhich.__new__.cache_clear()
        dviread.find_tex_file.cache_clear()
    except Exception:
        pass


_sanitize_tex_environment()

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.texmanager import TexManager
import subprocess
import warnings

def _apply_fallback_fonts():
    """Computer Modern mathtext without external LaTeX."""
    plt.rcParams.update({
        'text.usetex': False,
        'font.family': 'serif',
        'font.serif': ['DejaVu Serif', 'Times New Roman', 'serif'],
        'mathtext.fontset': 'cm',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'font.size': 14,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
        'axes.unicode_minus': False,
        'legend.framealpha': 0.9,
        'axes.axisbelow': True,
        'axes.grid': True,
        'grid.alpha': 0.3,
    })

def setup_latex_fonts():
    """
    Set up LaTeX fonts for Adobe Illustrator compatibility.
    
    Returns:
    --------
    bool: True if LaTeX is available and configured, False otherwise
    """
    
    # Reset to clean state
    plt.rcParams.update(plt.rcParamsDefault)
    
    # Set PDF backend
    if matplotlib.get_backend() != 'pdf':
        matplotlib.use('pdf')
    
    try:
        _sanitize_tex_environment()
        _reset_matplotlib_tex_caches()

        # Test LaTeX availability
        subprocess.check_output(['latex', '--version'], stderr=subprocess.DEVNULL)
        
        # Configure matplotlib for Adobe Illustrator compatibility
        plt.rcParams.update({
            # LaTeX configuration
            'text.usetex': True,
            'font.family': 'serif',
            'font.serif': ['Computer Modern Roman'],
            'font.sans-serif': ['Computer Modern Sans Serif'], 
            'font.monospace': ['Computer Modern Typewriter'],
            
            # Math fonts - Computer Modern only
            'mathtext.fontset': 'cm',
            'mathtext.rm': 'serif',
            
            # CRITICAL: LaTeX preamble without lmodern
            'text.latex.preamble': r'''
                \usepackage{amsmath}
                \usepackage{amsfonts}
                \usepackage{amssymb}
            ''',
            
            # Font embedding for Adobe compatibility
            'pdf.fonttype': 42,      # Type 1 fonts (Adobe compatible)
            'ps.fonttype': 42,       # Type 1 fonts for PostScript
            'svg.fonttype': 'none',  # Don't convert fonts to paths in SVG
            
            # Quality settings
            'font.size': 14,
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
            
            # Avoid font issues
            'axes.unicode_minus': False,
            
            # Visual improvements
            'legend.framealpha': 0.9,
            'axes.axisbelow': True,
            'axes.grid': True,
            'grid.alpha': 0.3,
        })

        # Verify matplotlib can actually render with usetex (not just latex --version)
        TexManager().get_text_width_height_descent(r'$\Delta V$', 14, False)
        return True
        
    except Exception as e:
        warnings.warn(f"LaTeX rendering unavailable: {e}. Using Computer Modern mathtext fallback.")
        _apply_fallback_fonts()
        return False

CMU_FONT_DIR = "/usr/share/texlive/texmf-dist/fonts/opentype/public/cm-unicode"
CMU_FONT_FILES = {
    "roman": "cmunrm.otf",
    "bold": "cmunbx.otf",
    "italic": "cmunti.otf",
    "bolditalic": "cmunbi.otf",
}
_CMU_FONTS_REGISTERED = False


def setup_illustrator_safe_fonts():
    """
    Configure genuine Computer Modern (CMU-Unicode OpenType) fonts WITHOUT
    invoking the external LaTeX/dvi pipeline.

    Why: matplotlib's `text.usetex=True` route renders true LaTeX Computer
    Modern glyphs, but those get embedded as Type 1 fonts with a custom
    ("Builtin") TeX encoding and no ToUnicode map. Adobe Illustrator often
    can't treat that text as live/editable text (font substitution, broken
    copy-paste, "missing font" warnings), even though it displays fine.

    This function instead registers the real CMU-Unicode OpenType Computer
    Modern fonts (shipped with TeX Live, `cm-unicode` package) directly with
    matplotlib's font manager and uses them via a custom mathtext fontset.
    With `pdf.fonttype = 42`, this embeds fully-subset CID/CFF OpenType
    fonts with a complete ToUnicode map -- verified with `pdffonts`/
    `pdftotext` to be selectable, searchable, and Illustrator-safe, with no
    external LaTeX subprocess required at all.

    Returns:
    --------
    bool: True if the CMU fonts were found and registered, False if it fell
    back to the generic fallback fonts.
    """
    global _CMU_FONTS_REGISTERED

    plt.rcParams.update(plt.rcParamsDefault)

    try:
        if not _CMU_FONTS_REGISTERED:
            for style_file in CMU_FONT_FILES.values():
                path = os.path.join(CMU_FONT_DIR, style_file)
                if not os.path.exists(path):
                    raise FileNotFoundError(f"CMU font not found: {path}")
                font_manager.fontManager.addfont(path)
            _CMU_FONTS_REGISTERED = True

        plt.rcParams.update({
            'text.usetex': False,
            'font.family': 'CMU Serif',
            'mathtext.fontset': 'custom',
            'mathtext.rm': 'CMU Serif',
            'mathtext.it': 'CMU Serif:italic',
            'mathtext.bf': 'CMU Serif:bold',

            # Full CID/CFF OpenType embedding with ToUnicode -- confirmed
            # Illustrator-safe (live, selectable text), unlike Type 3.
            'pdf.fonttype': 42,
            'ps.fonttype': 42,
            'svg.fonttype': 'none',

            'font.size': 14,
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
            'axes.unicode_minus': False,
            'legend.framealpha': 0.9,
            'axes.axisbelow': True,
            'axes.grid': True,
            'grid.alpha': 0.3,
        })
        return True
    except Exception as e:
        warnings.warn(f"CMU Illustrator-safe fonts unavailable: {e}. Using generic fallback.")
        _apply_fallback_fonts()
        return False


def latex_safe(latex_text, fallback_text=None, use_latex=True):
    """
    Safely format text for LaTeX or fallback rendering.
    
    Parameters:
    -----------
    latex_text : str
        Text with LaTeX formatting
    fallback_text : str, optional  
        Unicode fallback text
    use_latex : bool
        Whether LaTeX is available
        
    Returns:
    --------
    str: Properly formatted text
    """
    if use_latex:
        return latex_text
    else:
        if fallback_text:
            return fallback_text
        
        # Convert common LaTeX symbols to Unicode
        conversions = {
            r'\tau': 'τ', r'\xi': 'ξ', r'\varepsilon': 'ε', r'\epsilon': 'ε',
            r'\times': '×', r'\alpha': 'α', r'\beta': 'β', r'\gamma': 'γ',
            r'\delta': 'δ', r'\lambda': 'λ', r'\mu': 'μ', r'\nu': 'ν',
            r'\pi': 'π', r'\rho': 'ρ', r'\sigma': 'σ', r'\phi': 'φ',
            r'\omega': 'ω', r'\theta': 'θ', r'\kappa': 'κ', r'\eta': 'η',
            r'\mathrm{': '', '}': '', '$': '', r'\_': '_', r'\^': '^',
        }
        
        result = latex_text
        for latex_cmd, unicode_char in conversions.items():
            result = result.replace(latex_cmd, unicode_char)
        
        return result

# Set up style when module is imported
_USE_LATEX = setup_latex_fonts()

def get_latex_status():
    """Get whether LaTeX is available."""
    return _USE_LATEX

# Convenience functions for common formatting
def format_coupling(value):
    """Format coupling strength value for display."""
    if value == 0:
        return latex_safe(r'$\varepsilon = 0$', 'ε = 0', _USE_LATEX)
    elif value >= 1e-3:
        return latex_safe(f'$\\varepsilon = {value:.0e}$', f'ε = {value:.0e}', _USE_LATEX)
    else:
        return latex_safe(f'$\\varepsilon = {value:.0e}$', f'ε = {value:.0e}', _USE_LATEX)

def format_coupling_axis_label():
    """Format coupling strength axis label with units."""
    return latex_safe(r'$\varepsilon$ ($\times 10^{-4}$ a.u.)', 'ε (×10⁻⁴ a.u.)', _USE_LATEX)

def format_time(label='t'):
    """Format time variable."""
    if label == 't':
        return latex_safe(r'$t_{\mathrm{w}}$ (ps)', 't_w (ps)', _USE_LATEX)
    else:
        return latex_safe(f'${label}$ (ps)', f'{label} (ps)', _USE_LATEX)

def format_relaxation_time():
    """Format normalized relaxation time symbol."""
    return latex_safe(r'$\tilde{\tau}_{\mathrm{s}}$', 'τ̃_s', _USE_LATEX)

def format_fkt():
    """Format F(k,t) symbol."""
    return latex_safe(r'$F(k,t)$', 'F(k,t)', _USE_LATEX)

def format_material_time():
    """Format material time symbol.""" 
    return latex_safe(r'$\xi$ (ps)', 'ξ (ps)', _USE_LATEX)

if __name__ == "__main__":
    print("LaTeX Configuration Status:")
    print(f"LaTeX available: {_USE_LATEX}")
    print(f"Backend: {matplotlib.get_backend()}")
    print(f"Font type setting: {plt.rcParams['pdf.fonttype']}")
    print(f"Using LaTeX: {plt.rcParams['text.usetex']}")
