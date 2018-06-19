extensions [ palette table nw ]

; interactions = number of moves or imitations occurred after the execution of the Axelrod-Schelling model
; layoutsMap = map from layout strings to anonymous functions used to build them
globals [ interactions layoutsMap editDistanceMap standardColors colorMapping emptyCounter codeCounter codeTable
         plotInteractionCounter clusteringCoefficient density degree averagePathLength diameter ]

undirected-link-breed [ undirected-edges undirected-edge ]

turtles-own
[
  emptySite?        ;; if true, the site is not inhabited
  code              ;; site's cultural code
  omega             ;; site's average cultural overlap, USED FOR STATS PURPOSES ONLY. The main protocol model calls every time the appropriate function to get updated fresh values.
]

to setup
  clear-all
  resetRngSeed

  ; builds a map table from layoutChosen string to anonymous functions that build networks;
  ; this code is used to avoid nested ifelse
  set layoutsMap table:make
  table:put layoutsMap "spatially clustered network" [ -> setup-spatially-clustered-network ]
  table:put layoutsMap "preferential attachment" [ -> setupPreferentialAttachment ]
  table:put layoutsMap "Erdős–Rényi random network model" [ -> setupErdosRenyi ]
  table:put layoutsMap "Watts-Strogatz small world" [ -> setupWattsStrogatz ]
  table:put layoutsMap "Kleinberg model" [ -> setupKleinberg]

  ; builds a map from a chosen edit distance mechanism and the anonymous function used to compute it
  set editDistanceMap table:make
  table:put editDistanceMap "modified Hamming distance" [ [ a ] -> getModifiedHammingDistance a ]
  table:put editDistanceMap "cosine similarity" [ [ a ] -> getCosineDistance a ]

  set codeTable table:make

  set-default-shape turtles "circle"
  ; this commands do not affect subsequent random events
  with-local-randomness
  [
    ; list of standard colors to use if we want to color different values for a single cultural trait
    set standardColors sort sentence n-values 13 [i -> 13 + 10 * i] n-values 13 [i -> 15 + 10 * i]
    set colorMapping n-of q_value standardColors
  ]

  set interactions 0
  ;sets values used for stats purposes
  set codeCounter 0
  set clusteringCoefficient 0
  set density 0
  set degree 0
  set averagePathLength 0
  set diameter 0

  ; builds the chosen network layout
  let setupNetwork table:get layoutsMap layoutChosen
  run setupNetwork

  initializeNetwork

  with-local-randomness
  [
    set clusteringCoefficient mean [ nw:clustering-coefficient ] of turtles
    calculateAveragePathLength
    calculateDiameter
    set degree mean [count my-links] of turtles
    calculateDensity
    registerCodeStats
  ]

  set emptyCounter count turtles with [ emptySite? ]

  redoColor
  reset-ticks
end

to go
  clear-output
  with-local-randomness [ registerCodeStats ]
  axelrodSchelling
  ;redoColor
  with-local-randomness [ redoColor ]

  ; show interactions

  with-local-randomness [ collectAverageOverlap ]
  if interactions = 0
    [
      set plotInteractionCounter interactions
      stop
    ]

  ; we need another variable to save the old counter value, otherwise the plot will always show zero because interactions gets resetted at the end of the cycle
  set plotInteractionCounter interactions
  set interactions 0

  tick
end

to resetRngSeed
  ifelse fixedRandomSeed = true
  [
    ; to create replicable results use specific seeds
    ; random-seed 95199254
    random-seed customSeed
  ]
  [
    let seed new-seed
    print word "The current seed is: " seed
    random-seed seed
  ]
end

; Sets at least one node empty, other nodes are empty with a specific probability.
; Then non empty sites get a random cultural code.
to initializeNetwork
  ; at least one node must be empty
  ask turtle 0
  [
    set emptySite? true
    set omega 0
  ]

  ask turtles with [ who != 0 ]
  [
    ifelse random-float 1 <= emptyProbability
    [
      set emptySite? true
    ]
    [
      set emptySite? false
      ; chose traits at random uniformly
      set code n-values f_value [ i -> random q_value]
    ]
    set omega 0
  ]

  set emptyCounter count turtles with [ emptySite? ]
end

; this code comes from the Virus on a Network example in the Netlogo Library.
; As the name suggests, it builds spatially clustered networks
to setup-spatially-clustered-network

  create-turtles numberOfNodes
  [
    ; for visual reasons, we don't put any nodes *too* close to the edges
    setxy (random-xcor * 0.95) (random-ycor * 0.95)
  ]

  let num-links (averageNodeDegree * numberOfNodes) / 2
  while [count links < num-links ]
  [
    ask one-of turtles
    [
      let choice (min-one-of (other turtles with [not link-neighbor? myself])
                   [distance myself])
      if choice != nobody [ create-link-with choice ]
    ]
  ]
  ; make the network look a little prettier
  repeat 10
  [
    layout-spring turtles links 0.3 (world-width / (sqrt numberOfNodes)) 1
  ]
end

; this is the main Axelrod-Schelling model
to axelrodSchelling
  ask turtles with [ not emptySite?   ]
  [
    ;; if all neighbors are empty, directly move to another empty site
    ifelse not any? link-neighbors with [ not emptySite? ]
      [ move who ]
    [
      let peer one-of link-neighbors with [ not emptySite? ]

      ;; computes the kronecker's delta of each couple of items of the two corresponding item of the cultural code lists
      ;; then folds it by summing all the values
      let culturalOverlap getCulturalOverlap self peer

      ;; with probability equal to cultural overlap copies one trait of the selected peer
      ifelse random-float 1 <= culturalOverlap
      [
        if culturalOverlap != 1 [set interactions interactions + 1]

        let index random f_value
        let trait item index [ code ] of peer
        set code replace-item index code trait
      ]
      [
        ;; if the averageCulturalOverlap is lower than the threshold, move to an empty site
        let averageCulturalOverlap getAverageCulturalOverlap self
        if averageCulturalOverlap < T_threshold [ move who ]
      ]
    ]

  ]
end

; computes and returns the cultural overlap between n1 and n2
to-report getCulturalOverlap [ n1 n2 ]
  ;; computes the kronecker's delta of each couple of items of the two corresponding item of the cultural code lists
  ;; then folds it by summing all the values
  report (reduce + (map [ [a b] -> ifelse-value (a = b) [1] [0] ] [ code ] of n1 [ code ] of n2 )) / f_value
end

;; returns the average cultural overlap of n1 over its neighbors
to-report getAverageCulturalOverlap [ n1 ]
  let average 0
  ask link-neighbors with [ not emptySite? ]
  [
    set average average + getCulturalOverlap n1 self
  ]
  report average / count link-neighbors with [ not emptySite? ]
end

; move the turtle identified by turtleID in another random empty site
to move [ turtleID ]
  set interactions interactions + 1
  let newSite one-of turtles with [ emptySite? ]
  ask newSite [ set code [ code ] of turtle turtleID ]
  ;set [ code ] of newSite [ code ] of turtle turtleID
  ask newSite [ set emptySite? false ]
  ask turtle turtleID [set emptySite? true]
end

; computes the Euclidean norm of a list of values
to-report normalize [ n ]
  report sqrt ( reduce + ( map [ [ a ]  -> a * a ] n ) )
end

; computes the dot product between two lists of values
to-report dotProduct [ n1 n2 ]
  report ( reduce + ( map [ [ a b] -> a * b ] n1 n1) )
end

; computes the cosine similarity
to-report getCosineSimilarity [ n1 n2 ]
  report ( dotProduct n1 n2 ) / ( normalize n1 * normalize n2 + 0.000000001) ;; bias to prevent division by zero
end

to-report getCosineDistance [ codeList ]
  report getCosineSimilarity codeList n-values f_value [0]
end

;; Reports a modified version of the Hamming distance between the cultural code provided as argument and the list made up of f_value zeros.
;; Instead of counting the number of differente values between the code given and a string of zeros, it sums the different values in the correspondinc positions
;; and the positions themselves. We need this last correction to distinguish between codes like [0 0 9 9] and [9 9 0 0].
to-report getModifiedHammingDistance [ codeList ]
  let positions n-values f_value [ i -> i]
  report (reduce + (map [ [a b c] -> ifelse-value (a != b) [a + c] [0] ] codeList n-values f_value [0] positions))
end

;; colors the sites following these rules
;; if the site is empty, set the color to white;
;; otherwise, if the user does not choos to color code on a single trait basis, choose an appropriate color in the reverse HSV gradient range from lime green to red.
;; The more the cultural code is similar to the code made up of zeros, the "greener" is the site.
;; First and last rgb value picked from https://www.colorhexa.com/32cd32-to-ff0000
;;
;; The distance used to measure similarity between cultural codes is the one chosen by the user in the interface
;;
;; If instead, the user decides to color code values on a single trait basis, use randomly assigned colors to this specific trait values
to redoColor

  ifelse colorSingleTrait = false or q_value > 26
  [
    ask turtles
    [
      ifelse (emptySite?)
      [
        set color white
      ]
      [
        let distanceChosen table:get editDistanceMap editDistance
        let codeDistance  ( runresult  distanceChosen code )

        ;; minimum value possible is 0, maximum value is f_value times (q_value - 1) + (f_value - 1) (because of the positions)
        if editDistance = "modified Hamming distance"
        [ set color palette:scale-gradient [[50 205 50] [250 0 0]] codeDistance 0 (f_value * (q_value - 1) + (f_value - 1)) ]

        if editDistance = "cosine similarity"
        [ set color palette:scale-gradient [[50 205 50] [250 0 0]] codeDistance -1 1 ]
      ]
    ]
  ]
  [
    output-print (word "Trait chosen: " traitChosen ", out of " q_value " values")
    output-print "These are the color values assigned (for detailes, see Tools->Color Watches)"
    let i 0
    foreach colorMapping
    [
      c -> output-print (word "Color for value " i ": " c)
      set i i + 1
    ]

    ask turtles
    [
      ifelse (emptySite?)
      [ set color white ]
      [
        set color ( item (item traitChosen code) colorMapping )
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; This next section is a slightly modified version of the Preferential Attachment netlogo example present in the library;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; used for creating a new node
to make-node [old-node]
  create-turtles 1
  [
    if old-node != nobody
      [ create-link-with old-node
        ;; position the new node near its partner
        move-to old-node
        fd 8
      ]
  ]
end

;; This code is the heart of the "preferential attachment" mechanism, and acts like
;; a lottery where each node gets a ticket for every connection it already has.
;; While the basic idea is the same as in the Lottery Example (in the Code Examples
;; section of the Models Library), things are made simpler here by the fact that we
;; can just use the links as if they were the "tickets": we first pick a random link,
;; and than we pick one of the two ends of that link.
to-report find-partner
  report [one-of both-ends] of one-of links
end

to layout
  ;; the number 3 here is arbitrary; more repetitions slows down the
  ;; model, but too few gives poor layouts
  repeat 3 [
    ;; the more turtles we have to fit into the same amount of space,
    ;; the smaller the inputs to layout-spring we'll need to use
    let factor sqrt count turtles
    ;; numbers here are arbitrarily chosen for pleasing appearance
    layout-spring turtles links (1 / factor) (7 / factor) (1 / factor)
    display  ;; for smooth animation
  ]
  ;; don't bump the edges of the world
  let x-offset max [xcor] of turtles + min [xcor] of turtles
  let y-offset max [ycor] of turtles + min [ycor] of turtles
  ;; big jumps look funny, so only adjust a little each time
  set x-offset limit-magnitude x-offset 0.1
  set y-offset limit-magnitude y-offset 0.1
  ask turtles [ setxy (xcor - x-offset / 2) (ycor - y-offset / 2) ]
end

;; the UI and the code now use layout, this is for testing purposes when the previous procedure does not yield satisfactory results
to springLayout
  let factor sqrt count turtles
    if factor = 0 [ set factor 1 ]
    layout-spring turtles links (1 / factor) (14 / factor) (1.5 / factor)
    display
end

to-report limit-magnitude [number limit]
  if number > limit [ report limit ]
  if number < (- limit) [ report (- limit) ]
  report number
end

to setupPreferentialAttachment
  make-node nobody
  make-node turtle 0

  repeat numberOfNodes - 2
  [
    make-node find-partner
    layout
  ]

  repeat 20 [ layout ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; generates a network following the Erdős–Rényi model
to setupErdosRenyi
  nw:generate-random turtles undirected-edges numberOfNodes n-probability
  layout-circle sort turtles max-pxcor * 0.9
end

; generates a network following the Watts-Strogatz model
to setupWattsStrogatz
  nw:generate-watts-strogatz turtles undirected-edges numberOfNodes 2 n-probability
  ; William Thomas Tutte's layout
  layout-circle sort turtles max-pxcor * 0.9
  layout-tutte max-n-of (count turtles * 0.5) turtles [ count my-links ] links 12
end

;; generates a network following the Kleinberg model
to setupKleinberg
  let dim 0

  ; to simplify things, we round up the square root of selected numberOfNodes (the algorithm works by first building a lattice rows * columns)
  set dim round sqrt numberOfNodes
  set numberOfNodes dim * dim

  nw:generate-small-world turtles undirected-edges dim dim 2.0 false
  ; William Thomas Tutte's layout
  layout-circle sort turtles max-pxcor * 0.9
  layout-tutte max-n-of (count turtles * 0.5) turtles [ count my-links ] links 12
end

;; counts how many different cultural codes there are in the network. Moreover it stores this codes
;; in a table with their presence counter.
to registerCodeStats
  set codeCounter 0

  table:clear codeTable

  let inhabitedSites turtles with [ not emptySite? ]

  ask inhabitedSites
  [
    let currentCode code

    let equalsSites turtles with [ code = currentCode ]

    ; set codeCounter codeCounter + count equalsSites

    table:put codeTable currentCode count equalsSites

    set inhabitedSites turtles with [ code != currentCode ]
  ]

  set codeCounter table:length codeTable
end

;; computes and stores the average cultural overlap for each site
to collectAverageOverlap
  ask turtles with [not emptySite?]
  [
    if any? link-neighbors with [ not emptySite? ]
    [
      let averageCulturalOverlap getAverageCulturalOverlap self
      set omega averageCulturalOverlap
    ]
  ]
end

to colorOmegaLessThanOne
  ask turtles with [ emptySite? ]
  [
    set color white
  ]

  ask turtles with [ not emptySite?  and  omega = 1]
  [
    set color blue
  ]

  ask turtles with [ not emptySite? and omega < 1]
  [
    set color yellow
  ]
end

;; calculates the edge density
to calculateDensity
  set density 2 * (count links) / ( (count turtles) * (-1 + count turtles))
end

to calculateAveragePathLength
  let meanPathLength nw:mean-path-length

  ifelse meanPathlength = false
  [set averagePathLength "infinity"]
  [set averagePathLength meanPathLength]
end

to calculateDiameter
  ifelse member? false reduce [ [a b] -> sentence a b] [[ nw:distance-to myself ] of other turtles ] of turtles
  [set diameter "infinity"]
  [set diameter max [ max [ nw:distance-to myself ] of other turtles ] of turtles]
end

;; color communities found using the Louvain algorithm
to colorCommunities
  with-local-randomness
  [
    let communities nw:louvain-communities
    let colors sublist (sentence standardColors white grey) 0 (length communities)
    (foreach communities colors [ [community col] ->
      ask community [ set color col ] ])
  ]
end

;; color connected components
to colorComponents
  with-local-randomness
  [
    let components nw:weak-component-clusters
    let colors sublist (sentence shuffle standardColors white grey) 0 (length components)
    (foreach components colors [ [component col] ->
      ask component [ set color col ] ])
  ]
end

;; assigns a color to each cultural code and then color the turtles with that specific code
;; useful at the end of the simulation, when cultural codes will not be more than the cardinality of base-colors
to mapColorPopulation
  clear-output
  output-print "These are the color values assigned for the following codes (for color names, see Tools->Color Watches)"
  let currentCodes table:keys codeTable
  let colors sublist ( shuffle remove black (remove white base-colors) ) 0 (length currentCodes)
  (foreach currentCodes colors [ [curCode col] ->
    ask turtles with [code = curCode and not emptySite?] [ set color col ]
    output-print (word curCode ": " col)])
end
@#$#@#$#@
GRAPHICS-WINDOW
385
10
974
600
-1
-1
8.94
1
10
1
1
1
0
0
0
1
-32
32
-32
32
0
0
1
ticks
30.0

SLIDER
12
371
184
404
f_value
f_value
1
100
3.0
1
1
NIL
HORIZONTAL

SLIDER
194
371
366
404
q_value
q_value
1
100
3.0
1
1
NIL
HORIZONTAL

SLIDER
196
419
368
452
T_threshold
T_threshold
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
17
218
336
251
numberOfNodes
numberOfNodes
1
5000
400.0
1
1
NIL
HORIZONTAL

SLIDER
16
267
336
300
averageNodeDegree
averageNodeDegree
1
numberOfNodes - 1
6.0
1
1
NIL
HORIZONTAL

SLIDER
14
420
189
453
emptyProbability
emptyProbability
0
1
0.3
0.01
1
NIL
HORIZONTAL

BUTTON
35
73
116
106
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
151
73
241
106
go-once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
259
73
322
106
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

CHOOSER
14
11
319
56
layoutChosen
layoutChosen
"spatially clustered network" "preferential attachment" "Erdős–Rényi random network model" "Watts-Strogatz small world" "Kleinberg model"
0

TEXTBOX
23
308
323
374
averageNodeDegree is taken in consideration only for the spatially clustered layout, not for other networks
12
0.0
1

CHOOSER
15
480
263
525
editDistance
editDistance
"modified Hamming distance" "cosine similarity"
0

BUTTON
32
115
137
148
redo color
redoColor
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SWITCH
172
117
341
150
colorSingleTrait
colorSingleTrait
1
1
-1000

SLIDER
17
577
235
610
traitChosen
traitChosen
0
q_value - 1
0.0
1
1
NIL
HORIZONTAL

TEXTBOX
20
533
286
569
traitChosen used only if colorSingelTrait is On. In this case, editDistance is ignored
12
0.0
1

OUTPUT
1038
19
1719
153
14

SWITCH
17
627
240
660
fixedRandomSeed
fixedRandomSeed
1
1
-1000

PLOT
994
166
1412
434
Network status
time
# counters
0.0
100.0
0.0
100.0
true
true
"" ""
PENS
"interactions" 1.0 0 -16777216 true "" "plot plotInteractionCounter"
"number of different codes" 1.0 0 -2674135 true "" "plot codeCounter"

PLOT
1450
166
1885
434
Cultural codes
codes
number of nodes
0.0
30.0
0.0
100.0
true
true
"\n" "clear-plot\nlet keys table:keys codeTable\nwith-local-randomness\n[\n   let plotColors sentence n-values 13 [i -> 13 + 10 * i] n-values 13 [i -> 15 + 10 * i]\n   \n   (foreach keys [ k -> \n                     create-temporary-plot-pen (word k \"\")\n                     set-plot-pen-color one-of plotColors\n                     set-plot-pen-mode 1\n                     plotxy (position k keys + 1) table:get codeTable k ])\n]"
PENS

PLOT
992
465
1441
728
Cultural average overlaps
average overlap counter
counter
0.0
100.0
0.0
100.0
true
true
"" "clear-plot\nwith-local-randomness\n[\n   let plotColors sentence n-values 13 [i -> 13 + 10 * i] n-values 13 [i -> 15 + 10 * i]\n   let sites turtles with [ not emptySite? ]\n   let i 1\n   \n   ask sites\n   [\n      let tOmega omega\n      create-temporary-plot-pen (word tOmega \"\")\n      set-plot-pen-color one-of plotColors\n      set-plot-pen-mode 1\n      plotxy i count sites with [ omega = tOmega ]\n      set i i + 1.5\n      set sites turtles with [ omega != tOmega ]\n   ]\n]"
PENS

PLOT
1452
468
1885
730
Cultural average overlaps scatter plot
site ID
average overlap
0.0
100.0
0.0
1.5
true
false
"" "clear-plot\nwith-local-randomness\n[\n    ask turtles with [not emptySite?]\n    [\n       create-temporary-plot-pen \"scatter\"\n       set-plot-pen-color black\n       set-plot-pen-mode 2\n       plotxy who omega\n    ]\n]"
PENS

PLOT
438
614
918
803
ratio of sites with average overlap < 1
time
(#overlap<1) / #sites
0.0
100.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -2674135 true "" "plot count turtles with [not emptySite? and omega < 1] / count turtles with [not emptySite?]"

BUTTON
9
682
244
715
color sites with average overlap < 1
colorOmegaLessThanOne
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

TEXTBOX
9
663
301
681
shows in yellow sites with average omega < 1
12
0.0
1

MONITOR
992
737
1152
782
Clustering coefficient
clusteringCoefficient
17
1
11

MONITOR
1172
737
1330
782
Density
density
17
1
11

MONITOR
1345
738
1518
783
Average degree
degree
17
1
11

MONITOR
1530
737
1695
782
Average path length
averagePathLength
17
1
11

MONITOR
1713
737
1871
782
Diameter
diameter
17
1
11

BUTTON
33
158
137
191
redo-layout
springLayout
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
9
729
232
762
Color Louvain communities
colorCommunities
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

BUTTON
244
728
411
761
Color components
colorComponents
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
345
151
378
301
n-probability
n-probability
0
1
0.07
0.01
1
NIL
VERTICAL

TEXTBOX
158
159
338
229
n-probability used as a rewire probability for Watts-Strogatz and connection probability for Erdős–Rényi.
11
0.0
1

INPUTBOX
244
602
415
662
customSeed
0.0
1
0
Number

BUTTON
1730
114
1888
147
map codes to colors
with-local-randomness\n[\n   mapColorPopulation\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
1736
19
1886
109
this button is useful at the end of the simulation; it maps base-colors to cultural codes and colors the turtles
12
0.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
