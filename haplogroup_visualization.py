# -*- coding: utf-8 -*-
"""

haplogroup_visualization.py


Description:
This program creates a map displaying markers for different sample groups at different locations.
Clicking on a marker opens a sunburst chart showing the haplogroup frequencies present in that sample group.
Sunburst chart parts can be clicked to further zoom into subgroup frequencies of the haplogroups.
Each sunburst plot is labelled with the countries and age ranges the samples are from as well as the number of individuals included in that group.
The map also allows the user to switch between mt and Y chromosome data and filter for markers with individuals from certain age ranges.
The program opens the map automatically in the browser.


User-defined functions: create_Sunburst, creates interactive sunburst plot
Non-standard modules: numpy, pandas, plotly, folium, branca, and sklearn


Procedure:
    1. Define function for creating sunburst plots
    2. Input validation
        2.1 Check number of command-line arguments
        2.2 Verify input file name
    3. Data preprocessing
        3.1 Filter and clean dataset (remove irrelevant columns, correct country names, handle missing values)
        3.2 Add age coloumn for years in BC/AD
        3.3 Create separate tables for mt and Y data
        3.4 Group samples by chronological ranges
        3.5 Extract haplogroup subcategories
    4. Initialize map
    5. Apply K-Means clustering to group geographic locations
    6. Loop thorugh clusters to add markers and corresponding popup sunburst charts
    7. Save and open the interactive map in the browser


Input: AADR Annotations 2025.xlsx

Usage:
python haplogroup_visualization.py --cluster_Y (optional flag followed by an int) --cluster_mt (optional flag followed by an int) "AADR Annotations 2025.xlsx"


Version: 1.00
Date: 2025-03-25
Author: Lea Rachel Rieskamp
Git repository: https://github.com/RieskampLR/Population_Genetics_Project-2025.git


"""


# Library imports:

import os
import sys
import numpy as np
import pandas as pd
import argparse
import plotly.io as pio
import plotly.express as px
import folium
from folium.plugins import GroupedLayerControl, TagFilterButton
import branca
from sklearn.cluster import KMeans



###############################################################################

# Function for sunburst charts:
def create_Sunburst(subdataset):
    '''

    Parameters
    ----------
    subdataset : pandas dataframe
        Contains the data for the sunburst chart construction.

    Returns
    -------
    popup : folium.Popup object
        Contains an embedded sunburst chart in an HTML iframe.

    '''

    # Defining custom colours to ensure Haplogroups maintain same colour across charts
    color_map = {
    'A': '#4863A0', 'B': 'orange', 'C': '#FBBBB9', 'D': '#CC7A8B', 'E': '#FBE7A1',
    'G': '#BDF516', 'H': '#348781', 'I': '#654321', 'J': '#C83F49', 'K': '#9F000F',
    'L': '#004225', 'M': '#FEF250', 'N': '#46C7C7', 'O': '#550A35', 'P': 'indigo',
    'Q': '#C2E5D3', 'R': '#3B3131', 'T': '#667C26', 'U': '#9E7BFF', 'W': '#736F6E',
    'X': '#C83F49'
    }
    
    # Create sunburst
    sunburst = px.sunburst(subdataset, 
                           path=['first_letter', 'first_two_letters', 'first_three_letters', 'first_five_letters'],
                           maxdepth=3,
                           color='first_letter',
                           color_discrete_map=color_map)
    
    # Set hover-over labels
    sunburst.update_traces(hovertemplate="%{value}<br>%{percentRoot:.0%}")
    
    # Shift position to the right
    sunburst.update_layout(margin=dict(t=0, l=110, r=0, b=0))
   
    # Add legend
    popup_text = f'{cluster_countries}<br>{total_indivs} individuals<br>{cluster_range}'

    # Create sunburst html 
    html = (
        pio.to_html(sunburst, full_html=False, include_plotlyjs="cdn", config={'displaylogo': False,'displayModeBar': False}) 
        + f"<div style='position: absolute; left: 0px; top: 0px; font-family: Arial; font-size: 17px'>"
        + f"<h4>{popup_text}</h4></div>"
        )

    # Create iframe and popup and embed sunburst html
    iframe = branca.element.IFrame(html=html, width=380, height=260)
    popup = folium.Popup(iframe, max_width=380)
    
    return popup



###############################################################################

# Import data and error checks:

# Check of input argument number
if len(sys.argv) > 6:
    print("Error: Too many arguments. Please only provide the AADR excel annotations file and optionally requested cluster numbers.\nProgram terminated.")
    exit()

# Check of input file title
if sys.argv[-1] == "AADR Annotations 2025.xlsx":
    annotations = pd.read_excel(sys.argv[-1])
else:
    print('Error: Please provide the correct file ("AADR Annotations 2025.xlsx").\nProgram terminated')
    exit()



###################################################################################

# Cluster flags parser and flag input checks:

# Set up argument parser
parser = argparse.ArgumentParser(description="Specify cluster numbers for Y and mt")

# Add input file as a positional argument
parser.add_argument('input_file', type=str, help="AADR Annotations 2025.xlsx")

# Give error if no integer is given after a flag
try:
    # Define cluster arguments with optional flags
    parser.add_argument('--cluster_Y', type=int, default=150, help='Cluster number for Y (default: 150)')
    parser.add_argument('--cluster_mt', type=int, default=350, help='Cluster number for mt (default: 350)')
    # Parse arguments
    args = parser.parse_args()
except:
    print('Error: Please provide an integer between 5 and 500 following the --cluster flags.\nProgram terminated')
    exit()

# Check if cluster numbers are within restricted range
if not (5 <= args.cluster_Y <= 500):
    print("Error: Cluster number for Y must be between 5 and 500.\nProgram terminated.")
    exit()

if not (5 <= args.cluster_mt <= 500):
    print("Error: Cluster number for mt must be between 5 and 500.\nProgram terminated.")
    exit()



###############################################################################

# Modify annotation file data:
    
### Coloumns
# Filter out irrelevant coloumns
annotations = annotations.drop(columns=[
    "#",
    "Genetic ID",
    "Master ID",
    "Skeletal code",
    "Skeletal element",
    "Year data from this individual was first published [for a present-day individuals we give the data of the data reported here; missing GreenScience 2010 (Vi33.15, Vi33.26), Olalde2018 (I2657), RasmussenNature2010 (Australian)]",
    "Publication",
    "Method for Determining Date; unless otherwise specified, calibrations use 95.4% intervals from OxCal v4.4.2 Bronk Ramsey (2009); r5; Atmospheric data from Reimer et al (2020)", 
    "Date standard deviation in BP [OxCal sigma for a direct radiocarbon date, and standard deviation of the uniform distribution between the two bounds for a contextual date]",
    "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE",
    "Age at Death from physical anthropology",
    "Group ID",
    "Locality",
    "Pulldown Strategy",
    "Data source",
    "No. Libraries",
    "1240k coverage (taken from original pulldown where possible)",
    "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)",
    "SNPs hit on autosomal targets (Computed using easystats on HO snpset)",
    "Family ID and position within family",
    "Y haplogroup (manual curation in terminal mutation format)",
    "mtDNA coverage (merged data)",
    "mtDNA match to consensus if >2x (merged data)",
    "Damage rate in first nucleotide on sequences overlapping 1240k targets (merged data)",
    "Sex ratio [Y/(Y+X) counts] (merged data)",
    "Library type (minus=no.damage.correction, half=damage.retained.at.last.position, plus=damage.fully.corrected, ds=double.stranded.library.preparation, ss=single.stranded.library.preparation)",
    "Libraries",
    "ASSESSMENT",
    'ASSESSMENT WARNINGS (Xcontam interval is listed if lower bound is >0.005, "QUESTIONABLE" if lower bound is 0.01-0.02, "QUESTIONABLE_CRITICAL" or "FAIL" if lower bound is >0.02) (mtcontam confidence interval is listed if coverage >2 and upper bound is <0.'
    ])

# Rename coloumns
annotations = annotations.rename(columns={
    'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]':'Mean BP',
    'Political Entity':'Country',
    'Lat.':'Lat',
    'Long.':'Long',
    'Y haplogroup (manual curation in ISOGG format)': 'Y haplogroup',
    'mtDNA haplogroup if >2x or published':'mtDNA haplogroup'
    })


### Countries
# Change misspelled and long country names
annotations['Country'] = annotations['Country'].replace({'China ':'China', 
                                                         'Gernamy':'Germany', 
                                                         'Turkey ':'Turkey',
                                                         'Federated States of Micronesia':'Micronesia'
                                                         })


### BP
# Remove modern samples (BP=0)
annotations = annotations[annotations["Mean BP"] != 0]

# Add age coloumn (translation of BP to BC and AD based on Present=1950)
annotations["Age"] = np.where(
    annotations["Mean BP"] >= 1950, 
    (annotations["Mean BP"] - 1949).astype(str) + " BC",       # Convert to BC
    (1950 - annotations["Mean BP"]).astype(str) + " AD"        # Convert to AD
)


### Coordinates
# Remove missing Lat and Long coordinates
annotations = annotations[annotations["Lat"] != ".."]
annotations = annotations[annotations["Long"] != ".."]

# Remove NA values for Lat and Long
annotations = annotations.dropna(subset=["Lat", "Long"])

# Round Lat and Long values to 2 decimal points
annotations[["Lat", "Long"]] = annotations[["Lat", "Long"]].map(lambda x: round(x, 2))


### Haplogroups
# Convert haplogroup entries to strings
annotations["Y haplogroup"] = annotations["Y haplogroup"].astype(str)
annotations["mtDNA haplogroup"] = annotations["mtDNA haplogroup"].astype(str)

# Remove N/A and empty haplogroup entries
annotations = annotations[~annotations["Y haplogroup"].str.contains("n/a|na|NaN|not|Likely", na=False, regex=True)]
annotations = annotations[~annotations["mtDNA haplogroup"].str.contains("n/a|na|NaN|not|Likely", na=False, regex=True)]
annotations = annotations[annotations["Y haplogroup"].str.strip().ne('')]
annotations = annotations[annotations["mtDNA haplogroup"].str.strip().ne('')]

# Remove haplogroup naming additions
annotations["Y haplogroup"] = annotations["Y haplogroup"].str.replace(r"[+\/\(\)'~@\-or\s].*", "", regex=True)
annotations["mtDNA haplogroup"] = annotations["mtDNA haplogroup"].str.replace(r"[+\/\(\)'~@\-or\s].*", "", regex=True)


### Y and mt data and info
# Create separate tables for Y and mt data sets
annotations_Y = annotations.drop(columns=["mtDNA haplogroup"])
annotations_mt = annotations.drop(columns=["Y haplogroup"])

# Remove unknown haplogroups (storing the number of unknowns)
# Y
Y_unknowns_count = annotations_Y['Y haplogroup'].str.contains(r'\.\.').sum() 
annotations_Y = annotations_Y[annotations_Y["Y haplogroup"] != '..']
# mt
mt_unknowns_count = annotations_mt['mtDNA haplogroup'].str.contains(r'\.\.').sum()
annotations_mt = annotations_mt[annotations_mt["mtDNA haplogroup"] != '..']

# Most recent and oldest sample ages
# Y
youngest_Y = annotations_Y['Mean BP'].min()
oldest_Y = annotations_Y['Mean BP'].max()
# mt
youngest_mt = annotations_mt['Mean BP'].min()
oldest_mt = annotations_mt['Mean BP'].max()


### BP ranges
# Pre-grouping by BP range groups of 500
annotations['BP range'] = ((annotations['Mean BP'] // 500) * 500 + 1).astype(str) + '-' + ((annotations['Mean BP'] // 500 + 1) * 500).astype(str)
annotations_Y['BP range'] = ((annotations_Y['Mean BP'] // 500) * 500 + 1).astype(str) + '-' + ((annotations_Y['Mean BP'] // 500 + 1) * 500).astype(str)
annotations_mt['BP range'] = ((annotations_mt['Mean BP'] // 500) * 500 + 1).astype(str) + '-' + ((annotations_mt['Mean BP'] // 500 + 1) * 500).astype(str)

# Check for number of individuals per group
# Y
bp_counts = annotations_Y["BP range"].value_counts()
sorted_bp_counts = bp_counts.sort_index(key=lambda x: x.str.extract(r'(\d+)')[0].astype(int))
# mt
bp_counts = annotations_mt["BP range"].value_counts()
sorted_bp_counts = bp_counts.sort_index(key=lambda x: x.str.extract(r'(\d+)')[0].astype(int))

# Grouping by extended BP range groups
new_bp_ranges = [1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 11000, 44500]
bp_labels = ["1-1000", "1001-2000", "2001-3000", "3001-4000", "4001-5000", "5001-6000", "6001-7000", "7001-8000", "8001-11000", "11001-44500"]
annotations["BP range"] = pd.cut(annotations["Mean BP"], bins=new_bp_ranges, labels=bp_labels, right=True)
annotations_Y["BP range"] = pd.cut(annotations_Y["Mean BP"], bins=new_bp_ranges, labels=bp_labels, right=True)
annotations_mt["BP range"] = pd.cut(annotations_mt["Mean BP"], bins=new_bp_ranges, labels=bp_labels, right=True)

# Check for number of individuals per new group
# Y
bp_counts_Y = annotations_Y["BP range"].value_counts()
sorted_bp_counts_Y = bp_counts_Y.sort_index(key=lambda x: x.str.extract(r'(\d+)')[0].astype(int))
# mt
bp_counts_mt = annotations_mt["BP range"].value_counts()
sorted_bp_counts_mt = bp_counts_mt.sort_index(key=lambda x: x.str.extract(r'(\d+)')[0].astype(int))

# List of BP range categories for marker tags
bp_range_categories = annotations["BP range"].dropna().drop_duplicates().tolist()


### Haplogroup subcategories
# Subgrouping: 1, 2, 3, and 5 letters
# Y
annotations_Y['first_letter'] = annotations_Y['Y haplogroup'].str[:1]
annotations_Y['first_two_letters'] = annotations_Y['Y haplogroup'].str[:2]
annotations_Y['first_three_letters'] = annotations_Y['Y haplogroup'].str[:3]
annotations_Y['first_five_letters'] = annotations_Y['Y haplogroup'].str[:5]
# mt
annotations_mt['first_letter'] = annotations_mt['mtDNA haplogroup'].str[:1]
annotations_mt['first_two_letters'] = annotations_mt['mtDNA haplogroup'].str[:2]
annotations_mt['first_three_letters'] = annotations_mt['mtDNA haplogroup'].str[:3]
annotations_mt['first_five_letters'] = annotations_mt['mtDNA haplogroup'].str[:5]



###############################################################################

# Initialize map and set starting conditions:
m = folium.Map(location=[30,20], 
               tiles="Esri.WorldImagery", 
               zoom_start=2.5, 
               min_zoom=2, 
               max_zoom=7, 
               max_bounds=True)



###############################################################################

# Set group layers to separate Y and mt markers with interactive map buttons:
    
# Layers
fg_Y = folium.FeatureGroup(name='Y')
fg_mt = folium.FeatureGroup(name='mtDNA')
m.add_child(fg_Y)
m.add_child(fg_mt)

# Button
GroupedLayerControl(
    groups={'Haplogroups': [fg_Y, fg_mt]},
    collapsed=False,
).add_to(m)



###############################################################################

# K means clustering:

# Table of all coordinate combinations
lat_long_locations_Y = annotations_Y[[annotations_Y.columns[2], annotations_Y.columns[3]]].drop_duplicates()
lat_long_locations_mt = annotations_mt[[annotations_mt.columns[2], annotations_mt.columns[3]]].drop_duplicates()

# Convert lat_long to numpy arrays
coords_Y = np.array(lat_long_locations_Y[[annotations_Y.columns[2], annotations_Y.columns[3]]])
coords_mt = np.array(lat_long_locations_mt[[annotations_mt.columns[2], annotations_mt.columns[3]]])

# Set k (cluster numbers) based on user input
n_clusters_Y = args.cluster_Y
n_clusters_mt = args.cluster_mt

# Fit KMeans
kmeans_Y = KMeans(n_clusters=n_clusters_Y, random_state=42, n_init=10).fit(coords_Y)
kmeans_mt = KMeans(n_clusters=n_clusters_mt, random_state=42, n_init=10).fit(coords_mt)

# Add cluster labels
lat_long_locations_Y["Cluster"] = kmeans_Y.labels_
lat_long_locations_mt["Cluster"] = kmeans_mt.labels_



###############################################################################

# Add markers with popups to map - Y haplogroups:

# Loop through clusters
for cluster_id in range(n_clusters_Y):
    cluster_points = lat_long_locations_Y[lat_long_locations_Y["Cluster"] == cluster_id]

    # Get cluster center
    lat, long = kmeans_Y.cluster_centers_[cluster_id]

    # Subset cluster individuals
    subset = annotations_Y[
        annotations_Y[[annotations_Y.columns[2], annotations_Y.columns[3]]]
        .apply(tuple, axis=1)
        .isin(cluster_points[[annotations_Y.columns[2], annotations_Y.columns[3]]].apply(tuple, axis=1))
        ]
    total_indivs = subset.shape[0]

  # Find unique subset BP ranges
    bp_ranges = subset['BP range'].unique()
    bp_ranges = sorted(subset['BP range'].unique(), key=lambda x: int(x.split('-')[0]))
    cluster_first_bp = bp_ranges[0].split("-")[0]
    cluster_last_bp = bp_ranges[-1].split("-")[1]
    if cluster_first_bp == cluster_last_bp:
        cluster_range = f'{bp_ranges[0]} BP'
    else:
        cluster_range = f'{cluster_first_bp}-{cluster_last_bp} BP'
    countries = subset['Country'].unique()
    if len(countries) < 2:
        cluster_countries = countries[0]
    else:
        cluster_countries = f'{"<br>".join(countries)}'  
        
    # Create popup and add marker
    popup = create_Sunburst(subset)
    folium.CircleMarker(
        location=[lat, long],  
        radius=5,  
        color="black",
        fill=True,
        fill_opacity=1,
        fill_color="black",  
        popup=popup,
        tags=bp_ranges
    ).add_to(fg_Y)

    

###############################################################################

# Add markers with popups to map - mt haplogroups:

# Loop through clusters 
for cluster_id in range(n_clusters_mt):
    cluster_points = lat_long_locations_mt[lat_long_locations_mt["Cluster"] == cluster_id]

    # Get cluster center
    lat, long = kmeans_mt.cluster_centers_[cluster_id]

    # Subset cluster individuals
    subset = annotations_mt[
        annotations_mt[[annotations_mt.columns[2], annotations_mt.columns[3]]]
        .apply(tuple, axis=1)
        .isin(cluster_points[[annotations_mt.columns[2], annotations_mt.columns[3]]].apply(tuple, axis=1))
        ]
    total_indivs = subset.shape[0]
    
    # Find unique subset BP ranges
    bp_ranges = subset['BP range'].unique()
    bp_ranges = sorted(subset['BP range'].unique(), key=lambda x: int(x.split('-')[0]))
    cluster_first_bp = bp_ranges[0].split("-")[0]
    cluster_last_bp = bp_ranges[-1].split("-")[1]
    if cluster_first_bp == cluster_last_bp:
        cluster_range = f'{bp_ranges[0]} BP'
    else:
        cluster_range = f'{cluster_first_bp}-{cluster_last_bp} BP'
    countries = subset['Country'].unique()
    if len(countries) < 2:
        cluster_countries = countries[0]
    else:
        cluster_countries = f'{"<br>".join(countries)}'
        
    # Create popup and add marker
    popup = create_Sunburst(subset)
    folium.CircleMarker(
        location=[lat, long],  
        radius=5,  
        color="black",
        fill=True,
        fill_opacity=1,
        fill_color="black",  
        popup=popup,
        tags=bp_ranges
    ).add_to(fg_mt)



###############################################################################

# BP range filter:
    
# Sort BP range categories
bp_range_categories_sorted = sorted(bp_range_categories, key=lambda x: int(x.split('-')[0]))
# BP range tag-filter buttons
TagFilterButton(bp_range_categories_sorted).add_to(m)



###############################################################################

# Save and open map:

m.save("map.html")
os.system("map.html")


















