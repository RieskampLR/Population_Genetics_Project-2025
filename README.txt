README

A Novel Tool for visualising and exploring Haplogroup Frequencies in ancient Popula-tions on an interactive Map
Lea Rachel Rieskamp
25.03.2025


Overview
This project introduces an interactive tool for visualizing haplogroup distributions in ancient populations.
It allows users to explore Y and mtDNA haplogroup data across different time periods and geographic locations.
The tool presents data on a dynamic map with interactive sunburst charts displaying haplogroup frequencies.
It provides valuable insights into population dynamics, migration patterns, and haplogroup inheritance.


Methodology
The tool integrates geographical and chronological data and is built using Python (3.12.6).
Key libraries used:
- Folium: For creating the interactive map.
- Plotly: For generating sunburst charts.
- Branca: To embed HTMLs into Folium popups.
- Scikit-learn: For clustering samples geographically.
The tool supports:
Interactive map, navigable by zooming and dragging.
Popups displaying haplogroup frequencies and additional details.
Filtering by time ranges and switching between Y and mtDNA data.
Customizable cluster numbers via terminal flags.


Data and files
haplogroup_visualization.py: The Python script containing the code that generates the interactive map.
AADR_Annotations_2025.xlsx: The data file containing the Y and mtDNA data. Ensure this file is provided in its original form.


Installations & setup
Dependencies:
Ensure you have Python 3.12.6 and if necessary install the following libraries
- folium 0.19.5
- branca 0.8.1
- scikit-learn 1.6.1
- numpy 1.24.2
- pandas 2.0.1
- plotly 6.0.1
- argparse 1.4.0
Download the Data:
Download the AADR_Annotations_2025.xlsx file from the repository. Do not alter the file name or its contents.


Running the tool
Input: The input file is AADR_Annotations_2025.xlsx, which contains the necessary haplogroup data.
Output: The output is a dynamic web application that displays an interactive map of haplogroup distributions.
The map allows for navigation, filtering, and switching between Y and mtDNA data.
Terminal command: 
python haplogroup_visualization.py "AADR_Annotations_2025.xlsx"
This will launch an interactive map in your default web browser.


Customization of cluster numbers
If you wish to adjust the number of clusters for Y or mtDNA data, you can use the following flags when running the script:
python haplogroup_visualization.py --cluster_Y 200 --cluster_mt 400
This will set custom cluster numbers for the Y and mtDNA data (in the example above to 200 and 400).
The default cluster numbers are 150 (Y) and 350 (mt).
Both, either, or none of the flags can be applied.


Features
Navigable map: Zoom and drag to explore different regions.
Interactive popups: Click on a cluster to view additional information and a sunburst chart showing haplogroup frequencies.
Data filtering: Filter data by time ranges (e.g., 1-1000, 1001-2000, etc.).
Y/mtDNA toggle: Switch between Y and mtDNA data.
Sunburst chart interactions: Click to expand sunburst charts to see haplogroup frequencies in more detail.


Known issues & limitations
Performance:
Processing large datasets may take time.
Future improvements could include high-performance computing (HPC) support.
Compatibility:
The tool is currently designed to work with the provided AADR dataset,
and requires further modifications to support additional datasets.


Contact
For any questions, issues, or bug reports, please contact the author via mail:
lea.rieskamp@gmail.com



