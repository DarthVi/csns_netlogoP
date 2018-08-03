# CULTURAL DISSEMINATION AND SEGREGATION: AXELROD-SCHELLING MODEL

## WHAT IS IT?

This is a simulation about cultural dissemination and segregation, it implents the rules of the Axelrod-Schelling model, obtained by mixing Axelrod's model of dissemination of culture and Schelling's segregation model.

This project has been done as a project for the course of Complex Systems and Network Science at the Alma Mater University of Bologna and requires NetLogo >= 6.0.2.

## HOW IT WORKS
Nodes represent sites (houses, blocks of flats etc.), while edges represent neighbor relations. Each site can either be empty with probability _e_ (emptyProbability slider in the interface), or it can be inhabited by a single individual with probability (1 − _e_).
Each individual _i_ is associated with a non-empty site in the network and is characterized by a cultural code σ(_i_), which is a vector of length _F_ (`f_value` slider in the interface) where each element corresponds to a different cultural trait. 
Cultural traits are represented as integer values in the interval [0; _q_ − 1] (in the interface _q_ can be chosen by using the `q_value` slider). When a new individual is generated, each of the _F_ cultural traits in their cultural code σ(_i_) is picked at random uniformly within this interval. Once the network and the population have been generated, the model can be executed. 
At every simulation step, for each individual _i_ a random individual _j_
occupying one of the neighboring sites of _i_ is chosen. At this point, we compute
the cultural overlap ω<sub>_i_,_j_</sub> between _i_ and _j_ (see code for details about the formula). As the next step, _i_ randomly chooses a trait _k_ out of the _F_ cultural traits of _j_, and copies it with probability ω<sub>_i_,_j_</sub>.
If, on the other hand, imitation does not occur (with probability
(1 − ω<sub>_i_,_j_</sub>)), individual _i_ will compute the average cultural overlap with its entire neighborhood, which is defined as the average of the ω<sub>_i_,_j_</sub> values for each individual _j_ that is a neighbor of _i_. At this point if the average cultural overlap is lower than a fixed threshold _T_ (chosen by using the `T_threshold` slider in the interface), individual _i_ will decide to move to a random new site that is emtpy, leaving its previous site empty.

## HOW TO USE IT (interface explanations)

Select a network layout from the `layoutChosen` chooser then choose the number of nodes of the network from the `numberOfNodes` slider (it is advised to choose a value which is not too high, even if in theory you can go up to 1000; numerous experiments have been done with 400 nodes). If spatially clustered network is selected, you must decide the average node degree of the network, for other types of networks this settings will be ignored. For 2D lattices and the Kleinberg model, after choosing the numberOfNodes, the code will automatically adjust it to the closest appropriate value when the setup button is pressed, in order to be able to create a lattice of _k_ x _k_ nodes (this is done to simplify the coding section). 

Use the `n_probability` slider to define the rewire probability for the Watts-Strogatz layout or the connection probability for Erdős–Rényi network layout, as specified by the label. Then decide the values for the parameters of the Axelrod-Schelling model, which are the emptyProbability slider, f_value (cultural code length), q_value (each trait will be an integer in the interval between 0 and q_value - 1) and the T_threshold (tolerance threshold).

It is highly advised to activate the redo-layout button for the 2D Lattice, put the time slider on "faster" and deactivate the button when satisfied with the result.

After these steps, choose a coloring strategy (see the COLORING paragraph). You can also decide to use a fixed seed for the RNG, by putting fixedRandomSeed on On and writing the desired seed in the customSeed field. The availability of this option has been given in order to ensure replicability of the experiments.

Then press the setup button and the network will be generated. If you want, press redo-layout and deactivate it when desired. It is advised to use it with the 2D lattice, and the Erdős–Rényi layout. Spatially clustered network layout should produce an understandable and satisfiable layout on its own, without pressing redo-layout.
The Kleinberg model and the Watts-Strogatz model use William Thomas Tutte's layout, while
Preferential Attachment uses redo-layout in its own code by default.

Finally you can make the simulation run by pressing go (or go-once to run it for a single tick). At the end of the simulation it is advised to use the "map codes to colors" button (WARNING: this button may not work properly when the number of cultures over the population is higher than 27; in this case, use other methods). Other details about other button and plots will be explained in the following paragraphs.

See COLORING paragraph for the following buttons:

* Color Louvain communities
* Color components
* Color sites with average overlap < 1
* Maps codes to colors

### COLORING
Empty sites are always white, no matter the strategy chosen.
There are two main ways to color inhabited nodes while the simulation is running, in order to see how the model behaves:

1. compute the distance between the codes of the nodes we want to color and the code made up of f_value zeros; then assign accordingly a color in the reverse HSV gradient range from lime green (most similar to the string of zeros) to red (most different), see https://www.colorhexa.com/32cd32-to-ff0000 for getting information about the gradient scale used in the code.
To choose this strategy turn off the colorSingleTrait switch and pick the desired distance algorithm in the editDistance chooser (more info about these algorithms in the EDIT DISTANCE section).
WARNING: this strategy may produce results which are difficult to understand, mainly because similar codes have similar distances, which results in slightly different shades of the same color assigned to them, making them hard to distinguish. For this reason it is advised to use the "maps codes to colors" button at the end of the simulation.

2. turn on the colorSingleTrait switch and choose a trait to take in consideration from the traitChosen choose (ranging from 0 to f_value - 1). In this way the system will automatically assign random different colors to the possible q_value values of the chosen trait and color the nodes accordingly. The colors chosen by the algorithm will be shown in the output field. However keep in mind that the system will ignore this chose if q_value is higher than 26 (the number of colors available in order to make the nodes easily distinguishable) and will apply the first strategy. Moreover this strategy does not give you the opportunity to observe the overall network evolution, since only one trait is considered.

There are also additional buttons that influence coloring behaviours:

* redo-color: apply the coloring strategy chosen between the two previously discussed ones; it is usefule when you change idea about the strategy to use or about the trait to take in consideration
* Color Louvain communities: assign a random color (28 possible colors) to each Louvain community and color the nodes in the same community with the same assigned color; useful right after having pressed the setup button, before running the simulation, in order to understand how nodes are clustered togheter in hubs and clusters. The colors chosen by the algorithm will be shown in the output field along with the number of detected communities.
* Color components: colors each component with a different color (28 possible colors).
* color sites with average overlap < 1: colors with white the empty sites, blue for sites with average overlap equal to 1, yellow for sites with average overlap < 1; useful at the end of the simulation
* map codes to colors: empty sites will be shown in white. Assign other colors (27 possible colors) to each detected cultural codes and colors accordingly the nodes with the specific cultural codes. The colors chosen by the algorithm will be shown in the output field. It is advised to use this button at the end of the simulation. WARNING: this button may not work properly when the number of cultures over the population is higher than 27; in this case, use other methods.


### EDIT DISTANCE

When the colorSingleTrait switch is Off, coloring is based on an edit distance from the inhabited node's code and the code made up of f_value zeros (with f_value = 3, [0 0 0]), whereas empty sites are always white.
The main idea is that similar codes will have similar distance and will be colored accordingly.

The user can choose between two kinds of edit distance algorithms in the editDistance chooser: modified Hamming distance and cosine similarity.

* Modified Hamming distance: the original Hamming distance simply counts how many different values there are between two strings. This is not good enough for this simulation, because different codes like [1 1 0], [0 1 1] and [2 0 0] would have the same distance from zeros and be colored with the same code. For this reason the algorithm used in the code provided uses a modified version of it and makes a sum of the values that are different from zero, along with their position. For example [1 1 0] will have distance 3 (1 + 0 + 1 + 1, because the first 1 has position 0 and the second has position 1), [0 1 1] will have distance 5. Considering positions allows us to distinguish different codes that would have the same distance witouth this correction.

* Cosine similarity: it uses the cosine similarity algorithm with a bias ε of 0.000000001 at the denominator in order to prevent divisions by zero.

It is advised to used the modified Hamming distance and the "map codes to colors" button at the end of the simulation.

### OUTPUT FIELD

On the top right corner of the UI there is an output field that will be automaticalle filled with informations when needed. Here the seed of the RNG will be shown (in case fixedRandomSeed is Off) and also the number of detected cultures in the starting populations after pressing setup. Moreover if you hit the Color Louvain communities button, here will be shown the number of detected communities.
Other information that will be shown here: the mapping between colors and codes or traits whenever needed. However colors will be described with integers, since NegLogo's output-print command doesn't allow coloring text. Please refer to Tools -> Color Swatches to understand to what color a specific integer corresponds.

### PLOTS

ratio of sites with average overlap < 1: as the name suggests, it will show how the ratio of sites with average overlap < 1 (the number of these sites / all the inhabited sites) varies over time.

Network status: it shows how the number of interactions (imitations or nodes moving) varies and the number of different cultural codes in the population vary over time.

Cultural codes: it shows the histograms of detected cultural codes in the populations, with appropriate color coding in the legend; since at the beginning of the simulation there may be lots of cultural codes, the legend may not be fully displayable. The x-axis does not have a particular meaning, it only serves the purpose of separating the histograms.

Cultural average overlaps: it shows the different values cultural average overlap detected and draws their histograms, showing how many nodes hold each specific value.
The x-axis does not have a particular meaning, it only serves the purpose of separating the histograms.

Cultural average overlaps scatter plot: it shows a dot for each inhabited site with a specific cultural average overlap value.
X-axis: sites ID
Y-axis: cultural average overlap (the axis goes from 0 to 1.5 only for making the plot easy to read, the values themselves can reach a maximum of 1, they never go beyond that).

### MONITORS

On the bottom right side of the UI there are various monitors showing informations about the network: clustering coefficient, edge density, average node degree, average path length and network diameter.

## THINGS TO TRY

It is suggested to put the speed setting on "faster", especially when the number of chosen node is high.
Then follow the HOW TO USE IT steps; try for example the spatially clustered network layout with 150 nodes, q_value = 3, f_value = 3, then make different experiments with the empty probability and the threshold T. After getting accustomed to the UI with this values, try different configurations.

Most of the experiments have been done using the spatially clustered network layout, the 2D lattice and preferential attachment, because these layouts are the ones most similar to real cities layout (2D lattice for grid based cities and the spatially clustered for generic layouts) and online communities.

## NETLOGO FEATURES AND EXTENSIONS USED

Since NetLogo does not have a switch statement and ifelse can look ugly when too many of them are nested, the code makes use of the table extension. Each layout string is associated to an anonymous function which can be retrieved with a `table:get` using the layoutChosen value and then called.
Tables are also used to store the different cultural codes detected as keys of the table and their counter as content.

The palette estension is used to work more freely with colors. More specifically, it is used to define a color gradient from lime green to red via the function `palette:scale-gradient`.

The nw extension is used for building the Kleinberg layout, the Watts-Strogatz layout, the 2D lattice layout and the Erdős–Rényi layout (the other networks models and layout are built without using this extension). Moreover it is used to compute and show interesting data about the network, specifically the clustering coefficient, for detecting the Louvain communities, the average degree (in models different from the "spatially connected network" the user cannot choose the degre), the edge density, the average path length and the diameter. Some of these values (like the diameter) will be "infinity" if there isn't only one connected giant component.

## RELATED MODELS

See also models library -> Social Science -> Segregation.

The spatially clustered network layout code is the one used in models library -> Networks -> Virus on a Network

The preferential attachment layout code is a slightly modified version of the code in models library -> Networks -> Preferential Attachment

## CREDITS AND REFERENCES

Axelrod, Robert. “The dissemination of culture: A model with local convergence and global polarization.” Journal of conflict resolution 41.2 (1997):
203-226.

Schelling, Thomas C. “Dynamic models of segregation.” Journal of mathematical sociology 1.2 (1971): 143-186

Gracia-Lázaro, C., et al. “Residential segregation and cultural dissemination: An Axelrod-Schelling model.” Physical Review E 80.4 (2009): 046123

The spatially clustered network layout code is the one used in models library -> Networks -> Virus on a Network

The preferential attachment layout code is a slightly modified version of the code in models library -> Networks -> Preferential Attachment

## AUTHOR
Vito Vincenzo Covella a.k.a. "DarthVi"