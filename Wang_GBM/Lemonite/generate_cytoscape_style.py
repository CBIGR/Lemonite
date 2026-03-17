#!/usr/bin/env python3
"""
Generate a Cytoscape style (.xml) file for the LemonTree module-regulator network.

The style maps:
  - Node fill colour   → MegaGO cluster (discrete mapping on 'MegaGO_Cluster')
  - Node shape          → Node_Type  (ellipse for modules, diamond for regulators)
  - Node size           → Module_genes_count (passthrough / continuous for modules)
  - Node border colour  → Expression_significant (green = Yes, grey = No)
  - Node label          → shared name
  - Edge colour         → Category (Causal / Metabolic_pathway / Other)
  - Edge width          → Category
  - Edge target arrow   → ArrowShape column (DELTA / DOT / none)

Usage
-----
    python generate_cytoscape_style.py [--output FILE] [--clusters C1 C2 ...]

The resulting .xml can be imported into Cytoscape via
  File → Import → Styles from File ...

Requires only the Python standard library.
"""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
from pathlib import Path


# ── Colour palettes (must match the notebook) ────────────────────────────────
CLUSTER_PALETTE = [
    '#EE6677', '#4477AA', '#228833', '#AA3377', '#66CCEE',
    '#CCBB44', '#EE99AA', '#44BB99', '#BBCC33', '#AAAA00',
]
UNASSIGNED_COLOR = '#B0B0B0'
REGULATOR_COLOR  = '#FF8C00'

EDGE_COLORS = {
    'Causal':             '#2C3E50',
    'Metabolic_pathway':  '#E74C3C',
    'Other':              '#95A5A6',
}

EDGE_WIDTHS = {
    'Causal':            3.0,
    'Metabolic_pathway': 2.5,
    'Other':             1.0,
}

NODE_SHAPES = {
    'module':      'ELLIPSE',
    'TF':          'DIAMOND',
    'metabolite':  'DIAMOND',
    'lipid':       'HEXAGON',
}


# ── helpers ───────────────────────────────────────────────────────────────────

def _sub(parent: ET.Element, tag: str, text: str | None = None, **attribs) -> ET.Element:
    """Create a sub-element with optional text and attributes."""
    el = ET.SubElement(parent, tag, **attribs)
    if text is not None:
        el.text = str(text)
    return el


def _visual_property(parent, name, default=None, mapping_el=None):
    """Add a <visualProperty> element."""
    attribs = {'name': name}
    if default is not None:
        attribs['default'] = str(default)
    vp = ET.SubElement(parent, 'visualProperty', **attribs)
    if mapping_el is not None:
        vp.append(mapping_el)
    return vp


def _discrete_mapping(attr_name: str, attr_type: str, value_map: dict) -> ET.Element:
    """Build a <discreteMapping> element."""
    dm = ET.Element('discreteMapping', attributeName=attr_name, attributeType=attr_type.lower())
    for key, val in value_map.items():
        ET.SubElement(dm, 'discreteMappingEntry', attributeValue=str(key), value=str(val))
    return dm


def _passthrough_mapping(attr_name: str, attr_type: str) -> ET.Element:
    """Build a <passthroughMapping> element."""
    return ET.Element('passthroughMapping', attributeName=attr_name, attributeType=attr_type.lower())


def _continuous_mapping(attr_name: str, attr_type: str, points: list[dict]) -> ET.Element:
    """Build a <continuousMapping> with boundary points."""
    cm = ET.Element('continuousMapping', attributeName=attr_name, attributeType=attr_type.lower())
    for pt in points:
        ET.SubElement(cm, 'continuousMappingPoint', **{k: str(v) for k, v in pt.items()})
    return cm


# ── main builder ──────────────────────────────────────────────────────────────

def build_style(cluster_names: list[str] | None = None, style_name: str = 'LemonTree_MegaGO') -> ET.Element:
    """Return the root <vizmap> XML element."""

    if cluster_names is None:
        cluster_names = [f'Cluster_{i}' for i in range(1, 11)]

    # cluster → colour mapping
    cluster_color_map = {
        cl: CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)]
        for i, cl in enumerate(cluster_names)
    }
    cluster_color_map['Unassigned'] = UNASSIGNED_COLOR

    # Top-level element
    vizmap = ET.Element('vizmap', id='LemonTree', documentVersion='3.0')

    vs = ET.SubElement(vizmap, 'visualStyle', name=style_name)

    # ── network properties ────────────────────────────────────────────────
    net = ET.SubElement(vs, 'network')
    _visual_property(net, 'NETWORK_BACKGROUND_PAINT', '#FFFFFF')
    _visual_property(net, 'NETWORK_TITLE', style_name)

    # ── NODE properties ───────────────────────────────────────────────────
    node = ET.SubElement(vs, 'node')

    # Lock width/height so NODE_SIZE controls both
    dep = ET.SubElement(node, 'dependency', name='nodeSizeLocked', value='true')

    # Label = shared name
    _visual_property(node, 'NODE_LABEL', '',
                     mapping_el=_passthrough_mapping('shared name', 'String'))

    # Label font
    _visual_property(node, 'NODE_LABEL_FONT_SIZE', '12')
    _visual_property(node, 'NODE_LABEL_COLOR', '#000000')

    # Size mapped to gene count (modules only; regulators keep default)
    _visual_property(node, 'NODE_SIZE', '50',
                     mapping_el=_continuous_mapping('Module_genes_count', 'Integer', [
                         {'atlowerVal': '30', 'atupperVal': '30',
                          'equalVal':   '30', 'greaterVal':'30',
                          'lesserVal':  '30', 'value': '0'},
                         {'atlowerVal': '50', 'atupperVal': '50',
                          'equalVal':   '50', 'greaterVal':'50',
                          'lesserVal':  '50', 'value': '10'},
                         {'atlowerVal': '80', 'atupperVal': '80',
                          'equalVal':   '80', 'greaterVal':'80',
                          'lesserVal':  '80', 'value': '50'},
                         {'atlowerVal': '120', 'atupperVal': '120',
                          'equalVal':   '120', 'greaterVal':'120',
                          'lesserVal':  '120', 'value': '200'},
                     ]))

    # Fill colour by MegaGO cluster
    # Include regulator types so they also get coloured via the same column
    full_color_map = dict(cluster_color_map)
    full_color_map['regulator_TF'] = REGULATOR_COLOR
    full_color_map['regulator_metabolite'] = REGULATOR_COLOR
    full_color_map['regulator_lipid'] = REGULATOR_COLOR
    _visual_property(node, 'NODE_FILL_COLOR', UNASSIGNED_COLOR,
                     mapping_el=_discrete_mapping('MegaGO_Cluster', 'String', full_color_map))

    # Shape by Node_Type
    _visual_property(node, 'NODE_SHAPE', 'ELLIPSE',
                     mapping_el=_discrete_mapping('Node_Type', 'String', NODE_SHAPES))

    # Border width
    _visual_property(node, 'NODE_BORDER_WIDTH', '2')

    # Border colour: green if expression-significant, light grey otherwise
    _visual_property(node, 'NODE_BORDER_PAINT', '#CCCCCC',
                     mapping_el=_discrete_mapping('Expression_significant', 'String', {
                         'Yes': '#228833',
                         'No':  '#CCCCCC',
                     }))

    # Transparency
    _visual_property(node, 'NODE_TRANSPARENCY', '220')

    # Tooltip
    _visual_property(node, 'NODE_TOOLTIP', '',
                     mapping_el=_passthrough_mapping('shared name', 'String'))

    # ── EDGE properties ───────────────────────────────────────────────────
    edge = ET.SubElement(vs, 'edge')

    # Colour by Category
    _visual_property(edge, 'EDGE_STROKE_UNSELECTED_PAINT', EDGE_COLORS['Other'],
                     mapping_el=_discrete_mapping('Category', 'String', EDGE_COLORS))

    _visual_property(edge, 'EDGE_UNSELECTED_PAINT', EDGE_COLORS['Other'],
                     mapping_el=_discrete_mapping('Category', 'String', EDGE_COLORS))

    # Width by Category
    _visual_property(edge, 'EDGE_WIDTH', '1.0',
                     mapping_el=_discrete_mapping('Category', 'String', EDGE_WIDTHS))

    # Target arrow shape from ArrowShape column
    _visual_property(edge, 'EDGE_TARGET_ARROW_SHAPE', 'NONE',
                     mapping_el=_discrete_mapping('ArrowShape', 'String', {
                         'DELTA': 'DELTA',
                         'DOT':   'CIRCLE',
                     }))

    # Arrow colour matches stroke colour
    _visual_property(edge, 'EDGE_TARGET_ARROW_UNSELECTED_PAINT', EDGE_COLORS['Other'],
                     mapping_el=_discrete_mapping('Category', 'String', EDGE_COLORS))

    # Curve style
    _visual_property(edge, 'EDGE_BEND', '')
    _visual_property(edge, 'EDGE_CURVED', 'true')

    # Transparency
    _visual_property(edge, 'EDGE_TRANSPARENCY', '180')

    return vizmap


def prettify(element: ET.Element) -> str:
    """Return pretty-printed XML string."""
    rough = ET.tostring(element, encoding='unicode', xml_declaration=False)
    dom = minidom.parseString(rough)
    return dom.toprettyxml(indent='  ', encoding=None)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output .xml path (default: <output_dir>/cytoscape_style.xml)')
    parser.add_argument('--clusters', nargs='+', default=None,
                        help='Cluster names (e.g. Cluster_1 Cluster_2 …). '
                             'If omitted, Cluster_1 … Cluster_10 are used.')
    parser.add_argument('--style-name', type=str, default='LemonTree_MegaGO',
                        help='Name of the style inside Cytoscape (default: LemonTree_MegaGO)')
    parser.add_argument('--output-dir', type=str,
                        default='/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/'
                                'LemonTree/noProteomics_percentile2_divide_by_sum/'
                                'Networks/megaGO_exploration',
                        help='Directory to write the style file into (used when --output is not set)')
    args = parser.parse_args()

    vizmap = build_style(cluster_names=args.clusters, style_name=args.style_name)
    xml_str = prettify(vizmap)

    if args.output:
        out_path = Path(args.output)
    else:
        out_path = Path(args.output_dir) / 'cytoscape_style.xml'

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(xml_str, encoding='utf-8')
    print(f'Cytoscape style written to: {out_path}')


if __name__ == '__main__':
    main()
