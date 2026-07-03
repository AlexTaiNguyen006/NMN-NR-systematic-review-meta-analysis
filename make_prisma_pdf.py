#!/usr/bin/env python3
"""Convert PRISMA 2020 checklist CSV directly to a clean landscape PDF."""
import csv
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.units import cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
)
from reportlab.lib.enums import TA_CENTER

# Read CSV
rows = []
with open('submission_PHN/prisma_2020_checklist.csv') as f:
    reader = csv.reader(f)
    for row in reader:
        rows.append(row)

# Page setup
output = 'submission_PHN/prisma_2020_checklist.pdf'
doc = SimpleDocTemplate(
    output,
    pagesize=landscape(letter),
    topMargin=1.2*cm,
    bottomMargin=1.2*cm,
    leftMargin=1.2*cm,
    rightMargin=1.2*cm,
)

# Styles
styles = getSampleStyleSheet()
title_style = ParagraphStyle(
    'Title2', parent=styles['Title'],
    fontName='Times-Bold', fontSize=14, alignment=TA_CENTER,
    spaceAfter=6,
)
subtitle_style = ParagraphStyle(
    'Subtitle2', parent=styles['Normal'],
    fontName='Times-Italic', fontSize=10, alignment=TA_CENTER,
    spaceAfter=12,
)
header_style = ParagraphStyle(
    'Header', parent=styles['Normal'],
    fontName='Times-Bold', fontSize=8, leading=10,
)
cell_style = ParagraphStyle(
    'Cell', parent=styles['Normal'],
    fontName='Times-Roman', fontSize=7.5, leading=9.5,
)
section_style = ParagraphStyle(
    'Section', parent=styles['Normal'],
    fontName='Times-Bold', fontSize=8, leading=10,
)

# Build story
story = []
story.append(Paragraph('PRISMA 2020 Checklist', title_style))
story.append(Paragraph(
    'Mapping the Evidence Gap Between NMN and NR for Metabolic Outcomes: '
    'A Systematic Review, Transitivity Assessment, and Indirect Comparison '
    'Meta-Analysis',
    subtitle_style
))
story.append(Spacer(1, 6))

# Column widths
col_widths = [2.8*cm, 1.2*cm, 14*cm, 6*cm]

# Build table data with Paragraphs for text wrapping
table_data = []
section_rows = set()

for i, row in enumerate(rows):
    is_section = (
        len(row) >= 3
        and row[0].isupper()
        and row[1] == ''
        and row[2] == ''
    )

    if i == 0:
        table_data.append([
            Paragraph(row[0], header_style),
            Paragraph(row[1], header_style),
            Paragraph(row[2], header_style),
            Paragraph(row[3] if len(row) > 3 else '', header_style),
        ])
    elif is_section:
        section_rows.add(i)
        table_data.append([
            Paragraph(row[0], section_style),
            Paragraph('', cell_style),
            Paragraph('', cell_style),
            Paragraph('', cell_style),
        ])
    else:
        table_data.append([
            Paragraph(row[0], cell_style),
            Paragraph(row[1], cell_style),
            Paragraph(row[2] if len(row) > 2 else '', cell_style),
            Paragraph(row[3] if len(row) > 3 else '', cell_style),
        ])

# Create table
t = Table(table_data, colWidths=col_widths, repeatRows=1)

# Style
style_cmds = [
    ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
    ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#D9E2F3')),
    ('VALIGN', (0, 0), (-1, -1), 'TOP'),
    ('TOPPADDING', (0, 0), (-1, -1), 2),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
    ('LEFTPADDING', (0, 0), (-1, -1), 4),
    ('RIGHTPADDING', (0, 0), (-1, -1), 4),
]

# Shade section rows
for row_idx in section_rows:
    style_cmds.append(
        ('BACKGROUND', (0, row_idx), (-1, row_idx), colors.HexColor('#F2F2F2'))
    )
    style_cmds.append(
        ('SPAN', (0, row_idx), (-1, row_idx))
    )

t.setStyle(TableStyle(style_cmds))

story.append(t)
doc.build(story)
print(f'Created: {output}')
