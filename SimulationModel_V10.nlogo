extensions[gis csv]

breed[nodes node]
breed[birds bird]

nodes-own[node-id        ;;  each node has an id, which is from network analysis in R
          longitude      ;;  geogeraphic location
          latitude       ;;  geogeraphic location
          node-area       ;;  size of the water area
          node-country  ;;which country does the node belong to
          node-location  ;;  is the node wintering ground, passage or breeding ground?
          node-proportion-change ;;proportion that changes
          node-trend ;; is the node expanding or decreasing
          node-birds-aggregation-infected
          node-birds-aggregation
          node-birds-aggregation-immuned
          node-birds-aggregation-susceptible ; number of susceptible birds
          node-local-R0
          node-local-R0-contact
          node-local-R0-environment
          node-vrius-per-shed-rate
          invaded]

patches-own [patch-id                 ; corresponds to node-id
             patch-area               ; corresponds to node-area
             patch-country            ; corresponds to node-country
             patch-location           ; corresponds to node-location
             patch-proportion-change
             patch-trend
             patch-birds-aggregation ; number of birds
             patch-birds-aggregation-infected ; number of infected birds
             patch-birds-aggregation-immuned ; number of immunes birds
             patch-birds-aggregation-susceptible ; number of susceptible birds
             patch-density ; density of birds
             infection-probability ; probability for a susceptible birds becoming infected
             direct-infection-probability
             enivronment-infection-probability
             virus-in-environemnt-per-shed-rate ; amount of virus in shedding rate
             patch-migration-attractiveness    ; attractiveness
             patch-local-R0
             patch-local-R0-contact
             patch-local-R0-environment]

birds-own[node-current-location ; the node where the birds are staying
          node-target-location  ; the node where the birds are flying to
          patch-current-location; corresponds to node-current-location
          patch-target-location ; corresponds to node-target-location
          bird-location         ; corresponds to patch-location, indicating if the where is the brid
          target-link
          body-mass
          initial-body-mass
          day-not-flying

          ; infection variables
          susceptible?
          infected?
          immune?
          sick-time
          infection-route]

links-own[geo-distance             ;; distance between two nodes (km). it is calculate in R
          migration-resistance    ;; attributes from links. it is the un-likelihood for a bird to migrate over
          ;migration-pressure
          ]

globals [distance-farest-nodes-in-netlogo ;
         distance-farest-nodes-in-geo ; Geo-distance between two most north and south nodes, calculated in R (km)
         flying-distance-daily-in-netlogo ; converted from last line
         body-mass-consume-rate ; body-mass that was consumed per km?????? converted???
         %Infection ; proportion of infected birds in the population
         %Immune
         %CumulativeInfection  ; proportion of cumulative infected birds in the population
         infection-duration    ; days of infection
         ]
to setup
  clear-all
  import-network
  ask links [ ifelse Hide-Links? [hide-link][show-link]]
  ask birds [ ifelse Fly-Trail?   [pen-down][pen-up]]
  setup-birds-population
  setup-initial-infection
  update-prevalence
  setup-globals
  ;prepare-output-visited-links
  update-patches-variables
  update-display
  reset-ticks
end

to import-network
;; load shapefile so that I can convert geo coordinates to netlogo coordinates
  let nodes-geo-data gis:load-dataset "word-shape-files/nodes_gwfg.shp"
  ;; Set the world envelope to the union of all of our dataset's envelopes
  gis:set-world-envelope (gis:envelope-union-of (gis:envelope-of nodes-geo-data))
; choose netwok load | open files | read node metrics

; load nodes files
if Load-Migration-Network = "Network-Autumn-Complete" or Load-Migration-Network = "Network-Spring-Complete"
  [file-open "NodeList_NW_basic_Netlogo.txt"]
if Load-Migration-Network = "Network-Autumn-50-40Remove" or Load-Migration-Network = "Network-Spring-50-40Remove"
  [file-open "Nodelist_50_40Remove_Netlogo.txt"]
if Load-Migration-Network = "Network-Autumn-50-30Remove" or Load-Migration-Network = "Network-Spring-50-30Remove"
  [file-open "Nodelist_50_30Remove_Netlogo.txt"]
if Load-Migration-Network = "Network-Autumn-50-20Remove" or Load-Migration-Network = "Network-Spring-50-20Remove"
  [file-open "Nodelist_50_20Remove_Netlogo.txt"]
if Load-Migration-Network = "Network-Autumn-50-10Remove" or Load-Migration-Network = "Network-Spring-50-10Remove"
  [file-open "Nodelist_50_10Remove_Netlogo.txt"]
if Load-Migration-Network = "Network-Autumn-50-0Remove" or Load-Migration-Network = "Network-Spring-50-0Remove"
  [file-open "Nodelist_50_0Remove_Netlogo.txt"]
;read nodes files
while [not file-at-end?][
let items read-from-string (word "[" file-read-line "]")
create-nodes 1 [
set shape "circle"
set color green
set node-id    item 0 items
set latitude   item 2 items
set longitude  item 1 items
set node-area  item 12 items
set size       log node-area 10
set node-country item 9 items
set node-location   item 10 items
set node-proportion-change item 13 items
set node-trend item 14 items

if Node-Label = "Node-ID"  [set label node-id]
if Node-Label = "Node-Area"[set label node-area]
if Node-Label = "Node-Country"[set label node-country]
if Node-Label = "Node-Location"[set label node-location]
if Node-Label = "Node-Area-Trend"[set label node-trend]
if Node-Label = "Node-Area-Change"[set label node-proportion-change]

;; converting geo-distance to Netlogo-distance
let envelope gis:world-envelope ; [ xmin xmax ymin ymax ]
let xscale (max-pxcor - min-pxcor) / (item 1 envelope - item 0 envelope)
let yscale (max-pycor - min-pycor) / (item 3 envelope - item 2 envelope)
let netlogo-x (longitude - item 0 envelope) * xscale + min-pxcor
let netlogo-y (latitude - item 2 envelope) * yscale + min-pycor
setxy netlogo-x netlogo-y
;show list netlogo-x netlogo-y
                ]        ]
file-close

; load autumn link files
if Load-Migration-Network = "Network-Autumn-Complete"
[file-open "Edgelist_NW_basic_Netlogo.txt"   read-autumn-network]
if Load-Migration-Network = "Network-Autumn-50-40Remove"
[file-open "Edgelist_NW_50_40Remove_Netlogo.txt" read-autumn-network]
if Load-Migration-Network = "Network-Autumn-50-30Remove"
[file-open "Edgelist_NW_50_30Remove_Netlogo.txt" read-autumn-network]
if Load-Migration-Network = "Network-Autumn-50-20Remove"
[file-open "Edgelist_NW_50_20Remove_Netlogo.txt" read-autumn-network]
if Load-Migration-Network = "Network-Autumn-50-10Remove"
[file-open "Edgelist_NW_50_10Remove_Netlogo.txt" read-autumn-network]
if Load-Migration-Network = "Network-Autumn-50-0Remove"
[file-open "Edgelist_NW_50_0Remove_Netlogo.txt" read-autumn-network]
; load spring link files
 if Load-Migration-Network = "Network-Spring-Complete"
[file-open "Edgelist_NW_basic_Netlogo.txt"   read-spring-network]
if Load-Migration-Network = "Network-Spring-50-40Remove"
[file-open "Edgelist_NW_50_40Remove_Netlogo.txt" read-spring-network]
if Load-Migration-Network = "Network-Spring-50-30Remove"
[file-open "Edgelist_NW_50_30Remove_Netlogo.txt" read-spring-network]
if Load-Migration-Network = "Network-Spring-50-20Remove"
[file-open "Edgelist_NW_50_20Remove_Netlogo.txt" read-spring-network]
if Load-Migration-Network = "Network-Spring-50-10Remove"
[file-open "Edgelist_NW_50_10Remove_Netlogo.txt" read-spring-network]
if Load-Migration-Network = "Network-Spring-50-0Remove"
[file-open "Edgelist_NW_50_0Remove_Netlogo.txt" read-spring-network]
; set nodes metrics to patches
ask nodes [set patch-id                node-id
           set patch-area              node-area
           set patch-country           node-country
           set patch-location          node-location
           set patch-proportion-change node-proportion-change
           set patch-trend             node-trend
           set pcolor green
           set invaded           0
          ]
end
to read-autumn-network ;read autumn link files
 while [not file-at-end?][
 let items read-from-string (word "[" file-read-line "]")
 ask get-node (item 0 items)[create-link-to get-node (item 1 items) [
    set geo-distance item 2 items
      set migration-resistance  exp (read-from-string geo-distance * 0.1)
    ;set migration-resistance  read-from-string geo-distance ^ 2
    set thickness 0.1] ]]
 file-close
end

to read-spring-network ; read spring link files
while [not file-at-end?][
 let items read-from-string (word "[" file-read-line "]")
 ask get-node (item 1 items)[create-link-to get-node (item 0 items) [
    set geo-distance item 2 items
    set migration-resistance  exp (read-from-string geo-distance * 0.1)
    set thickness 0.1] ]]
 file-close
end
to-report get-node [id]
report one-of nodes with [node-id = id]
end

to setup-birds-population
create-birds Population-Size [
  set shape  "default"
  set size    1
  set color   yellow
  ;set initial-body-mass random-normal Mean-Body-Mass Body-Mass-SD
  set body-mass random-normal Mean-Body-Mass Body-Mass-SD
  set day-not-flying 0
  ;set body-mass initial-body-mass
  if Load-Migration-Network = "Network-Autumn-Complete" or Load-Migration-Network = "Network-Autumn-50-40Remove" or Load-Migration-Network = "Network-Autumn-50-30Remove" or
     Load-Migration-Network = "Network-Autumn-50-20Remove"  or Load-Migration-Network = "Network-Autumn-50-10Remove" or Load-Migration-Network = "Network-Autumn-50-0Remove"
    [set node-current-location one-of nodes with-max [pycor]
     set patch-current-location [patch-here] of node-current-location]

  if Load-Migration-Network = "Network-Spring-Complete" or Load-Migration-Network = "Network-Spring-50-40Remove" or Load-Migration-Network = "Network-Spring-50-30Remove" or
     Load-Migration-Network = "Network-Spring-50-20Remove"  or Load-Migration-Network = "Network-Spring-50-10Remove" or Load-Migration-Network = "Network-Spring-50-0Remove"
    [set node-current-location one-of nodes with-min [pycor]
     set patch-current-location [patch-here] of node-current-location]

     move-to patch-current-location
     set node-target-location  node-current-location
     set patch-target-location patch-current-location
     set bird-location patch-location
     ]
end

to setup-initial-infection
  ask birds [get-susceptible
            set sick-time 0]
  ask n-of (Initial-Infected-Ratio * Population-Size) birds [get-sick]
  ask nodes [set virus-in-environemnt-per-shed-rate Initial-Virus-in-Environment]
end

to get-susceptible
set susceptible? true
set infected?    false
set immune?      false
set sick-time    0
end

to get-sick
set susceptible? false
set infected?    true
set immune?      false
end

to get-immune
set susceptible? false
set infected?    false
set immune?      true
set sick-time    0
end

to update-prevalence
set %Infection (count birds with [infected?]  / count birds)
set %Immune    (count birds with [immune?] / count birds)
set %CumulativeInfection (%Infection + %Immune)
ask nodes with [invaded = 0]
[if node-local-R0 > 1 [set invaded 1]]
end

to setup-globals
;set Body-Mass-Accumulate-Rate  46.3 ; g/day from Jasper Madsen 2006
set distance-farest-nodes-in-netlogo  88.68432863925128 ;; observer> ask node 97 [show distance node 0]
set distance-farest-nodes-in-geo  4509.709 ;; km. calculate in R #GeoDistances[nrow(GeoDistances),1]
;set Flying-Speed  483.84  ;; km/d =5.6 m/s, from Andrea Kölzsch et al. 2016 same speed for autumn and spring due to avoiding complesity that caused by speed difference
set flying-distance-daily-in-netlogo (Flying-Speed / ( distance-farest-nodes-in-geo / distance-farest-nodes-in-netlogo ))
set body-mass-consume-rate  Body-Mass-Comsume-Rate * ( distance-farest-nodes-in-netlogo / distance-farest-nodes-in-geo)  ;;0.3g/km is calculated
;;  with consumption rate 8.9 kJ/km and energy content 27 kJ/g from Jesper Madsen et al. 2006
set infection-duration 7  ;; 7days infection duration is general
end

to update-patches-variables
ask nodes [set patch-birds-aggregation          count birds-here
           set node-birds-aggregation  patch-birds-aggregation
           set patch-birds-aggregation-infected count birds-here with [ infected? ]
           set node-birds-aggregation-infected patch-birds-aggregation-infected
           set patch-birds-aggregation-immuned  count birds-here with [ immune? ]
           set node-birds-aggregation-immuned  patch-birds-aggregation-immuned
           set patch-birds-aggregation-susceptible  count birds-here with [ susceptible? ]
           set node-birds-aggregation-susceptible patch-birds-aggregation-susceptible
           set patch-density                    patch-birds-aggregation / patch-area
           set patch-migration-attractiveness 1 / exp patch-density
           set node-vrius-per-shed-rate virus-in-environemnt-per-shed-rate
           ;set patch-migration-attractiveness  exp patch-density / (1 + exp patch-density)
           ;set patch-migration-attractiveness  (exp 6.9 + (patch-density * -0.14)) / (1 + (exp 6.9 + (patch-density * -0.14)))
           ifelse patch-birds-aggregation != 0
             [set infection-probability             beta * (patch-birds-aggregation-infected + virus-in-environemnt-per-shed-rate) / patch-birds-aggregation
              set direct-infection-probability      beta *   patch-birds-aggregation-infected                                       / patch-birds-aggregation
              set enivronment-infection-probability beta *                                     virus-in-environemnt-per-shed-rate  / patch-birds-aggregation
              if  infection-probability >= 1    [set infection-probability 0.99]]

              [set infection-probability            0
               set direct-infection-probability     0
               set enivronment-infection-probability 0]]
;; update R0 for nodes
ask nodes[
ifelse patch-birds-aggregation = 0
  [set patch-local-R0-contact      0
   set patch-local-R0-environment  0]
  [ifelse patch-birds-aggregation-infected = 0
    [set patch-local-R0-contact 0]
    [set patch-local-R0-contact beta * infection-duration * patch-birds-aggregation-susceptible / patch-birds-aggregation]
   ifelse virus-in-environemnt-per-shed-rate = 0
    [set patch-local-R0-environment 0]
    [set patch-local-R0-environment beta * infection-duration * ((1 / Virus-Decaying-Rate) - 1) * patch-birds-aggregation-susceptible / patch-birds-aggregation]
  ]

set patch-local-R0             patch-local-R0-contact + patch-local-R0-environment

set node-local-R0-contact     patch-local-R0-contact
set node-local-R0-environment patch-local-R0-environment
set node-local-R0             patch-local-R0]

end

to update-virus-in-environment
ask nodes [set virus-in-environemnt-per-shed-rate virus-in-environemnt-per-shed-rate * ( 1 - Virus-Decaying-Rate) + patch-birds-aggregation-infected * ( 1 - Virus-Decaying-Rate)]
end

to update-display
ask birds[set color ifelse-value infected? [red][ifelse-value immune? [grey][yellow]]]
ask nodes with [color = green]
[if patch-birds-aggregation-infected >= Habitat-Contamination-Threshold * Population-Size [set color red]]
;ask nodes with [color = green]
;[if patch-birds-aggregation >= 1 [set color orange]]
end


;;go;;
;;go;;
;;go;;
;;go;;
;;go;;
;;go;;
;;go;;
;;go;;
;;go;;
;;go;;
;;go;;


to go
  ifelse when-stop-simulation = "all-arrive" [
if Load-Migration-Network = "Network-Autumn-Complete" or Load-Migration-Network = "Network-Autumn-50-40Remove" or Load-Migration-Network = "Network-Autumn-50-30Remove" or
     Load-Migration-Network = "Network-Autumn-50-20Remove"  or Load-Migration-Network = "Network-Autumn-50-10Remove" or Load-Migration-Network = "Network-Autumn-50-0Remove"
   ;[if all? birds [not infected?]  [stop]]
  [if all? birds [ (bird-location  = "Wintering") and (patch-here = patch-target-location)]  [stop]]
if Load-Migration-Network = "Network-Spring-Complete" or Load-Migration-Network = "Network-Spring-50-40Remove" or Load-Migration-Network = "Network-Spring-50-30Remove" or
     Load-Migration-Network = "Network-Spring-50-20Remove"  or Load-Migration-Network = "Network-Spring-50-10Remove" or Load-Migration-Network = "Network-Spring-50-0Remove"
    [if all? birds [ (bird-location  = "Breeding") and (patch-here = patch-target-location)]   [stop]]

  ]

  [ifelse when-stop-simulation = "no-infected" [
    if Load-Migration-Network = "Network-Autumn-Complete" or Load-Migration-Network = "Network-Autumn-50-40Remove" or Load-Migration-Network = "Network-Autumn-50-30Remove" or
     Load-Migration-Network = "Network-Autumn-50-20Remove"  or Load-Migration-Network = "Network-Autumn-50-10Remove" or Load-Migration-Network = "Network-Autumn-50-0Remove"
   ;[if all? birds [not infected?]  [stop]]
    [if all? birds [ (bird-location  = "Wintering") and (patch-here = patch-target-location) and (not infected?)]  [stop]]
if Load-Migration-Network = "Network-Spring-Complete" or Load-Migration-Network = "Network-Spring-50-40Remove" or Load-Migration-Network = "Network-Spring-50-30Remove" or
     Load-Migration-Network = "Network-Spring-50-20Remove"  or Load-Migration-Network = "Network-Spring-50-10Remove" or Load-Migration-Network = "Network-Spring-50-0Remove"
    [if all? birds [ (bird-location  = "Breeding") and (patch-here = patch-target-location) and (not infected?)]   [stop]]
  ][]]

ask birds [ ifelse patch-here = patch-target-location
              [synchronize-location
               transmission-process
               accumulate-body-mass

               ;if body-mass >= initial-body-mass * ( 1 + Migration-Body-Mass-Threshold ) [
                if body-mass >= Mean-Body-Mass * ( 1 + Migration-Body-Mass-Threshold ) [
                looking-for-target
                move]]

            [if infected? [ set sick-time sick-time + 1
                            if sick-time >= infection-duration [get-immune]]
                            move]
             ]
ask links [ ifelse Hide-Links? [hide-link][show-link]]
ask birds [ ifelse Fly-Trail?   [pen-down][pen-up]]
update-display
update-prevalence
update-patches-variables
update-virus-in-environment
  if Save-Visited-Links? [update-output-visited-links]
  ;plot [patch-birds-aggregation] of one-of patches with [patch-id = "1"]
tick
end

to transmission-process
  ;; the local variable infection-prob is for reducing the computational time, then the infection-probability only need to be calculated once per tick
  if susceptible? [ let infection-prob             infection-probability
                    let direct-infection-prob      direct-infection-probability
                    let enivronment-infection-prob enivronment-infection-probability
                    let dice random-float 1
                    if dice <= direct-infection-prob [get-sick
                                                      set infection-route 1]
                    if (direct-infection-prob <= dice) and (dice <= infection-prob) [get-sick
                                                                                     set infection-route 2]
                    ]
  if infected? [set sick-time sick-time + 1
                if sick-time >= infection-duration [get-immune]]
end

to accumulate-body-mass
set body-mass ( body-mass + Body-Mass-Accumulate-Rate )
set day-not-flying day-not-flying + 1
end

to looking-for-target
  if any? [out-link-neighbors] of node-current-location [let target-links-all [my-out-links] of node-current-location
                                                         ifelse any? target-links-all with [read-from-string geo-distance >= (distance-farest-nodes-in-geo / (1 + Numer-of-Stops))][
                                                         let target-links target-links-all with [read-from-string geo-distance >= (distance-farest-nodes-in-geo / (1 + Numer-of-Stops))]
                                                         ;
                                                         ifelse random-float 1 <= (1 - random-migrate-probability)
                                                          [set target-link max-one-of target-links [migration-pressure]
                                                           set node-target-location  [end2] of target-link
                                                           set patch-target-location [patch-here] of node-target-location]
                                                         [ set node-target-location [end2] of one-of target-links
                                                           set patch-target-location [patch-here] of node-target-location]
                                                           ]

                                                           [let target-links target-links-all
                                                           ;
                                                         ifelse random-float 1 <= 1 - random-migrate-probability
                                                          [set target-link max-one-of target-links [migration-pressure]
                                                           set node-target-location  [end2] of target-link
                                                           set patch-target-location [patch-here] of node-target-location]
                                                         [ set node-target-location [end2] of one-of target-links
                                                           set patch-target-location [patch-here] of node-target-location]
                                                           ]
                                                           ]
end


to move
  face patch-target-location
  ifelse distance patch-target-location <= flying-distance-daily-in-netlogo
    [move-to patch-target-location
     set body-mass body-mass - ( (distance patch-target-location) * body-mass-consume-rate )]
    [fd flying-distance-daily-in-netlogo
     set body-mass body-mass - (flying-distance-daily-in-netlogo * body-mass-consume-rate )]
end

to synchronize-location
set bird-location patch-location
set node-current-location  node-target-location
set patch-current-location patch-target-location
end

to-report migration-pressure
report ([patch-migration-attractiveness] of end2 - [patch-migration-attractiveness] of end1 ) / migration-resistance
end


to update-output-visited-links
  if Load-Migration-Network = "Network-Autumn-Complete" [file-open "Output_Visited_links_Autumn-Complete.csv"]
  if Load-Migration-Network = "Network-Autumn-50-40Remove" [file-open "Output_Visited_links_Autumn-50-40Remove.csv"]
  if Load-Migration-Network = "Network-Autumn-50-30Remove" [file-open "Output_Visited_links_Autumn-50-30Remove.csv"]
  if Load-Migration-Network = "Network-Autumn-50-20Remove" [file-open "Output_Visited_links_Autumn-50-20Remove.csv"]
  if Load-Migration-Network = "Network-Autumn-50-10Remove" [file-open "Output_Visited_links_Autumn-50-10Remove.csv"]
  if Load-Migration-Network = "Network-Autumn-50-0Remove" [file-open "Output_Visited_links_Autumn-50-0Remove.csv"]

  if Load-Migration-Network = "Network-Spring-Complete" [file-open "Output_Visited_links_Spring-Complete.csv"]
  if Load-Migration-Network = "Network-Spring-50-40Remove" [file-open "Output_Visited_links_Spring-50-40Remove.csv"]
  if Load-Migration-Network = "Network-Spring-50-30Remove" [file-open "Output_Visited_links_Spring-50-30Remove.csv"]
  if Load-Migration-Network = "Network-Spring-50-20Remove" [file-open "Output_Visited_links_Spring-50-20Remove.csv"]
  if Load-Migration-Network = "Network-Spring-50-10Remove" [file-open "Output_Visited_links_Spring-50-10Remove.csv"]
  if Load-Migration-Network = "Network-Spring-50-0Remove" [file-open "Output_Visited_links_Spring-50-0Remove.csv"]

  ask birds with [node-current-location != node-target-location]
            [file-type (word (ticks + 1) ",")
              file-type (word self ",")
            file-type (word [node-id] of node-current-location ",")
            file-print (word [node-id] of node-target-location ",")]
file-close
end
@#$#@#$#@
GRAPHICS-WINDOW
439
10
1255
827
-1
-1
8.0
1
10
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
days
30.0

BUTTON
341
46
407
79
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
254
90
317
123
NIL
go
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
318
90
381
123
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
1

CHOOSER
3
90
253
135
Load-Migration-Network
Load-Migration-Network
"Network-Autumn-Complete" "Network-Autumn-50-40Remove" "Network-Autumn-50-30Remove" "Network-Autumn-50-20Remove" "Network-Autumn-50-10Remove" "Network-Autumn-50-0Remove" "Network-Spring-Complete" "Network-Spring-50-40Remove" "Network-Spring-50-30Remove" "Network-Spring-50-20Remove" "Network-Spring-50-10Remove" "Network-Spring-50-0Remove"
11

CHOOSER
2
44
173
89
Node-Label
Node-Label
"Node-ID" "Node-Area" "Node-Country" "Node-Location" "Node-Area-Trend" "Node-Area-Change"
0

SLIDER
176
150
348
183
Population-Size
Population-Size
0
15000
10000.0
10
1
NIL
HORIZONTAL

SLIDER
4
182
176
215
Mean-Body-Mass
Mean-Body-Mass
0
5000
2025.0
5
1
NIL
HORIZONTAL

SLIDER
176
182
348
215
Body-Mass-SD
Body-Mass-SD
0
500
150.0
5
1
NIL
HORIZONTAL

SLIDER
6
426
187
459
Initial-Infected-Ratio
Initial-Infected-Ratio
0
1
0.01
0.01
1
NIL
HORIZONTAL

SLIDER
6
459
192
492
Virus-Decaying-Rate
Virus-Decaying-Rate
0
1
0.11
0.01
1
NIL
HORIZONTAL

SLIDER
8
527
295
560
Habitat-Contamination-Threshold
Habitat-Contamination-Threshold
0
1
1.0E-4
0.001
1
NIL
HORIZONTAL

SLIDER
5
216
262
249
Migration-Body-Mass-Threshold
Migration-Body-Mass-Threshold
0
1
0.15
0.01
1
NIL
HORIZONTAL

SLIDER
6
393
178
426
beta
beta
0
2
0.09
0.01
1
NIL
HORIZONTAL

SLIDER
6
249
237
282
Flying-Speed
Flying-Speed
0
1000
486.0
1
1
NIL
HORIZONTAL

PLOT
9
689
178
809
Infection Dynamic
Days
Precentage
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"Infected" 1.0 0 -2674135 true "" "plot %Infection"
"Immuned" 1.0 0 -7500403 true "" "plot %Immune"
"Cumulative Infected" 1.0 0 -955883 true "" "plot %CumulativeInfection"

INPUTBOX
6
315
155
375
Body-Mass-Accumulate-Rate
38.0
1
0
Number

INPUTBOX
157
315
306
375
Body-Mass-Comsume-Rate
0.5
1
0
Number

SWITCH
1
10
130
43
Hide-Links?
Hide-Links?
1
1
-1000

SWITCH
131
10
257
43
Fly-Trail?
Fly-Trail?
1
1
-1000

SLIDER
7
493
236
526
Initial-Virus-in-Environment
Initial-Virus-in-Environment
0
1000
0.0
10
1
NIL
HORIZONTAL

SLIDER
6
281
235
314
random-migrate-probability
random-migrate-probability
0
0.05
0.001
0.0001
1
NIL
HORIZONTAL

SLIDER
4
149
176
182
Numer-of-Stops
Numer-of-Stops
0
20
7.0
1
1
NIL
HORIZONTAL

CHOOSER
175
45
338
90
when-stop-simulation
when-stop-simulation
"all-arrive" "no-infected" "no-restrict"
2

PLOT
180
689
354
810
direct-vs-environment
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"direct" 1.0 0 -7500403 true "" "plot count birds with [infection-route = 1]"
"environment" 1.0 0 -955883 true "" "plot count birds with [infection-route = 2]"
"two-routes" 1.0 0 -6459832 true "" "plot (count birds with [infection-route = 1] + count birds with [infection-route = 2])"

SWITCH
259
10
437
43
Save-Visited-Links?
Save-Visited-Links?
1
1
-1000

PLOT
8
563
179
686
Node_Local_R0
NIL
NIL
0.0
5.0
0.0
1.0
true
false
"" ""
PENS
"node1" 1.0 0 -16777216 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"1\"]"
"node3" 1.0 0 -7500403 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"3\"]"
"node 40" 1.0 0 -13345367 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"40\"]"
"node62" 1.0 0 -955883 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"62\"]"
"node 70" 1.0 0 -6459832 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"70\"]"
"node 86" 1.0 0 -1184463 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"86\"]"
"node 22" 1.0 0 -10899396 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"22\"]"
"node 33" 1.0 0 -13840069 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"33\"]"
"node 90" 1.0 0 -14835848 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"90\"]"
"node 98" 1.0 0 -2674135 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"98\"]"
"node 96" 1.0 0 -955883 true "" "plot [node-local-R0] of one-of nodes with [node-id = \"96\"]"

PLOT
181
564
355
686
Node_Local_R0_Direct
NIL
NIL
0.0
10.0
0.0
3.0
true
false
"" ""
PENS
"node1" 1.0 0 -16777216 true "" "plot [node-local-R0-contact] of one-of nodes with [node-id = \"1\"]"
"node3" 1.0 0 -7500403 true "" "plot [node-local-R0-contact] of one-of nodes with [node-id = \"3\"]"
"node40" 1.0 0 -2674135 true "" "plot [node-local-R0-contact] of one-of nodes with [node-id = \"40\"]"
"node62" 1.0 0 -955883 true "" "plot [node-local-R0-contact] of one-of nodes with [node-id = \"62\"]"
"node70" 1.0 0 -6459832 true "" "plot [node-local-R0-contact] of one-of nodes with [node-id = \"70\"]"
"node86" 1.0 0 -1184463 true "" "plot [node-local-R0-contact] of one-of nodes with [node-id = \"86\"]"
"node98" 1.0 0 -10899396 true "" "plot [node-local-R0-contact] of one-of nodes with [node-id = \"98\"]"

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
NetLogo 6.0.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="NW_Autumn_stop_no_restrict" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>count birds with [infected?]</metric>
    <metric>count birds with [susceptible?]</metric>
    <metric>count birds with [immune?]</metric>
    <metric>count birds with [infection-route = 1]</metric>
    <metric>count birds with [infection-route = 2]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "98"]</metric>
    <metric>mean               [day-not-flying] of birds</metric>
    <enumeratedValueSet variable="Body-Mass-SD">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Virus-in-Environment">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Flying-Speed">
      <value value="486"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Node-Label">
      <value value="&quot;Node-ID&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Load-Migration-Network">
      <value value="&quot;Network-Autumn-Complete&quot;"/>
      <value value="&quot;Network-Autumn-50-20Remove&quot;"/>
      <value value="&quot;Network-Autumn-50-0Remove&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Numer-of-Stops">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Infected-Ratio">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean-Body-Mass">
      <value value="2025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Accumulate-Rate">
      <value value="38"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat-Contamination-Threshold">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Migration-Body-Mass-Threshold">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Size">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="when-stop-simulation">
      <value value="&quot;no-restrict&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fly-Trail?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Comsume-Rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random-migrate-probability">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hide-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Save-Visited-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Virus-Decaying-Rate">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="NW_Autumn_stop_no_restrict" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>count birds with [infected?]</metric>
    <metric>count birds with [susceptible?]</metric>
    <metric>count birds with [immune?]</metric>
    <metric>count birds with [infection-route = 1]</metric>
    <metric>count birds with [infection-route = 2]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "98"]</metric>
    <metric>mean               [day-not-flying] of birds</metric>
    <enumeratedValueSet variable="Body-Mass-SD">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Virus-in-Environment">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Flying-Speed">
      <value value="486"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Node-Label">
      <value value="&quot;Node-ID&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Load-Migration-Network">
      <value value="&quot;Network-Autumn-Complete&quot;"/>
      <value value="&quot;Network-Autumn-50-40Remove&quot;"/>
      <value value="&quot;Network-Autumn-50-30Remove&quot;"/>
      <value value="&quot;Network-Autumn-50-20Remove&quot;"/>
      <value value="&quot;Network-Autumn-50-10Remove&quot;"/>
      <value value="&quot;Network-Autumn-50-0Remove&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Numer-of-Stops">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Infected-Ratio">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean-Body-Mass">
      <value value="2025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Accumulate-Rate">
      <value value="38"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat-Contamination-Threshold">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Migration-Body-Mass-Threshold">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Size">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="when-stop-simulation">
      <value value="&quot;no-restrict&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fly-Trail?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Comsume-Rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random-migrate-probability">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hide-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Save-Visited-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Virus-Decaying-Rate">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="SA_beta_NW_Autumn_stop_no_restrict" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>count birds with [infected?]</metric>
    <metric>count birds with [susceptible?]</metric>
    <metric>count birds with [immune?]</metric>
    <metric>count birds with [infection-route = 1]</metric>
    <metric>count birds with [infection-route = 2]</metric>
    <enumeratedValueSet variable="Body-Mass-SD">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Virus-in-Environment">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Flying-Speed">
      <value value="486"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.01"/>
      <value value="0.2"/>
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Node-Label">
      <value value="&quot;Node-ID&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Load-Migration-Network">
      <value value="&quot;Network-Autumn-Complete&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Numer-of-Stops">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Infected-Ratio">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean-Body-Mass">
      <value value="2025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Accumulate-Rate">
      <value value="38"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat-Contamination-Threshold">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Migration-Body-Mass-Threshold">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Size">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="when-stop-simulation">
      <value value="&quot;no-restrict&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fly-Trail?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Comsume-Rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random-migrate-probability">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hide-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Save-Visited-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Virus-Decaying-Rate">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="NW_Spring_stop_no_restrict" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>count birds with [infected?]</metric>
    <metric>count birds with [susceptible?]</metric>
    <metric>count birds with [immune?]</metric>
    <metric>count birds with [infection-route = 1]</metric>
    <metric>count birds with [infection-route = 2]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-birds-aggregation-infected] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-birds-aggregation] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-local-R0] of one-of nodes with [node-id = "98"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "1"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "2"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "3"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "4"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "5"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "6"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "7"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "8"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "9"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "10"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "11"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "12"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "13"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "14"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "15"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "16"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "17"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "18"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "19"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "20"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "21"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "22"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "23"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "24"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "25"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "26"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "27"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "28"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "29"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "30"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "31"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "32"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "33"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "34"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "35"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "36"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "37"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "38"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "39"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "40"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "41"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "42"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "43"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "44"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "45"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "46"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "47"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "48"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "49"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "50"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "51"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "52"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "53"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "54"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "55"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "56"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "57"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "58"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "59"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "60"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "61"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "62"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "63"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "64"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "65"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "66"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "67"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "68"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "69"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "70"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "71"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "72"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "73"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "74"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "75"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "76"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "77"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "78"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "79"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "80"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "81"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "82"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "83"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "84"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "85"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "86"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "87"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "88"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "89"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "90"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "91"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "92"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "93"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "94"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "95"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "96"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "97"]</metric>
    <metric>[node-local-R0-contact] of one-of nodes with [node-id = "98"]</metric>
    <metric>mean               [day-not-flying] of birds</metric>
    <enumeratedValueSet variable="Body-Mass-SD">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Virus-in-Environment">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Flying-Speed">
      <value value="486"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Node-Label">
      <value value="&quot;Node-ID&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Load-Migration-Network">
      <value value="&quot;Network-Spring-Complete&quot;"/>
      <value value="&quot;Network-Spring-50-20Remove&quot;"/>
      <value value="&quot;Network-Spring-50-0Remove&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Numer-of-Stops">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Infected-Ratio">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean-Body-Mass">
      <value value="2025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Accumulate-Rate">
      <value value="38"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat-Contamination-Threshold">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Migration-Body-Mass-Threshold">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Size">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="when-stop-simulation">
      <value value="&quot;no-restrict&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fly-Trail?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Comsume-Rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random-migrate-probability">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hide-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Save-Visited-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Virus-Decaying-Rate">
      <value value="0.1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="SA_beta_noEnv_NW_Autumn_stop_no_restrict" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>count birds with [infected?]</metric>
    <metric>count birds with [susceptible?]</metric>
    <metric>count birds with [immune?]</metric>
    <metric>count birds with [infection-route = 1]</metric>
    <metric>count birds with [infection-route = 2]</metric>
    <enumeratedValueSet variable="Body-Mass-SD">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Virus-in-Environment">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Flying-Speed">
      <value value="486"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.01"/>
      <value value="0.2"/>
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Node-Label">
      <value value="&quot;Node-ID&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Load-Migration-Network">
      <value value="&quot;Network-Autumn-Complete&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Numer-of-Stops">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Infected-Ratio">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean-Body-Mass">
      <value value="2025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Accumulate-Rate">
      <value value="38"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat-Contamination-Threshold">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Migration-Body-Mass-Threshold">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Size">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="when-stop-simulation">
      <value value="&quot;no-restrict&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fly-Trail?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Comsume-Rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random-migrate-probability">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hide-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Save-Visited-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Virus-Decaying-Rate">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="SA_eta_NW_Autumn_stop_no_restrict" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="150"/>
    <metric>count birds with [infected?]</metric>
    <metric>count birds with [susceptible?]</metric>
    <metric>count birds with [immune?]</metric>
    <metric>count birds with [infection-route = 1]</metric>
    <metric>count birds with [infection-route = 2]</metric>
    <enumeratedValueSet variable="Body-Mass-SD">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Virus-in-Environment">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Flying-Speed">
      <value value="486"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Node-Label">
      <value value="&quot;Node-ID&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Load-Migration-Network">
      <value value="&quot;Network-Autumn-Complete&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Numer-of-Stops">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Initial-Infected-Ratio">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Mean-Body-Mass">
      <value value="2025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Accumulate-Rate">
      <value value="38"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Habitat-Contamination-Threshold">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Migration-Body-Mass-Threshold">
      <value value="0.15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Population-Size">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="when-stop-simulation">
      <value value="&quot;no-restrict&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Fly-Trail?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Body-Mass-Comsume-Rate">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="random-migrate-probability">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Hide-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Save-Visited-Links?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Virus-Decaying-Rate">
      <value value="0.01"/>
      <value value="0.02"/>
      <value value="0.03"/>
      <value value="0.04"/>
      <value value="0.05"/>
      <value value="0.06"/>
      <value value="0.07"/>
      <value value="0.08"/>
      <value value="0.09"/>
      <value value="0.1"/>
      <value value="0.11"/>
      <value value="0.12"/>
      <value value="0.13"/>
      <value value="0.14"/>
      <value value="0.15"/>
      <value value="0.16"/>
      <value value="0.17"/>
      <value value="0.18"/>
      <value value="0.19"/>
      <value value="0.2"/>
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
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
