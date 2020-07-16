rm(list = ls())

{library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggnetwork)
library(network)
library(deSolve)
library(igraph)
library(GGally)}

select <- dplyr::select


# ---------------------------------------
#     Test
# ---------------------------------------

highschool

df <- data.frame(from = c(1, 1), to = c(2, 3))


network <- network(df, matrix.type = "edgelist")

ggraph(network) + 
  geom_edge_link() + 
  geom_node_point()


network <- network(highschool %>% slice(1:10), directed = TRUE, matrix.type = "edgelist")
# network <- graph_from_data_frame(highschool %>% slice(1:10), directed = TRUE)

ggraph(network) + 
  geom_edge_link() + 
  geom_node_point()




# ---------------------------------------
#     Plotting pars
# ---------------------------------------
colors <- tableau_color_pal(palette = "Tableau 10", 
                            type = c("regular",
                                     "ordered-sequential", 
                                     "ordered-diverging"), direction = 1)(5)[c(1, 3, 5)]

pie(rep(1, length(colors)), col = colors)


# ---------------------------------------
#     Loading data
# ---------------------------------------
ReactionHistory.long <- read_csv("../output/reactionHistory.csv")
InfectionReactions <- ReactionHistory.long %>% filter(ReactionType == "Infection") %>%
  mutate(Source = str_extract(ReactionNodes, "(?<=\\[)[0-9]+") %>% as.numeric(),
         Target = str_extract(ReactionNodes, "[0-9]+(?=\\])") %>% as.numeric()) %>%
  select(-ReactionNodes)

InfectionNbrs <- InfectionReactions %>% select(Source) %>%
  distinct() %>%
  left_join(
    InfectionReactions %>% group_by(Source) %>%
      summarize(InfectionsBySource = n())    
  ) %>%
  left_join(
    InfectionReactions %>% select(t, Target) %>%
      group_by(Target) %>% filter(t == min(t)) %>%
      rename(Source = Target,
             TimeOfInfection = t)
  ) %>%
  mutate(TimeOfInfection = replace_na(TimeOfInfection, 0)) %>%
  arrange(TimeOfInfection) %>% as.data.frame()

InfectionNbrs %>%
  ggplot() +
  geom_point(aes(x = TimeOfInfection, y = InfectionsBySource)) +
  ylim(c(0, 7))


SimSummary.long <- read_csv("../output/summarizedDynamicState.csv") %>%
  pivot_longer(-c(t, iter), names_to = "State", values_to = "Count")

SimSummary.long %>%
  ggplot() +
  geom_line(aes(x = t, y = Count, group = State, color = factor(State, levels = c("S", "I", "R")))) +
  scale_color_manual(name = "State", values = colors)


SimSummary.long %>%
  filter(t < 20, State != "S") %>%
  ggplot() +
  geom_line(aes(x = t, y = Count, group = State, color = factor(State, levels = c("S", "I", "R")))) +
  scale_color_manual(name = "State", values = colors) +
  scale_y_log10()


# Growth rate
lm.fit <- lm(log(Count) ~ t, data = SimSummary.long %>% filter(t < 18, State == "I"))
R0 <- 2.5
gamma <- 0.2
(R0 - 1) * gamma

1 + lm.fit$coefficients[2] / gamma

InfectionNbrs %>%
  filter(TimeOfInfection < 20) %>%
  pull(InfectionsBySource) %>%
  summary()



Sim <- read_csv("../output/dynamicState.csv")
IC <- Sim %>% slice(1) %>% select(-1) %>%
  pivot_longer(cols = everything(), names_to = "Node")




network.tib <- read_csv("../output/network.csv") %>%
  mutate(Node2 = as.numeric(Node2)) %>%
  arrange(Node1, Node2)

node.names <- network.tib %>%
  pivot_longer(cols = everything()) %>%
  filter(!is.na(value)) %>% select(value) %>%
  distinct() %>% arrange(value) %>%
  pull(value)

network <- network.initialize(n = length(node.names), directed = FALSE)

for (l in 1:nrow(network.tib)){
  node1 <- network.tib %>% slice(l) %>% pull(Node1)
  node2 <- network.tib %>% slice(l) %>% pull(Node2)
  
  if ( !is.na(node2)){
    add.edge(network, node1, node2)
  }
}

network %v% "degree" <- degree(network)
# network %v% "state" <- Sim %>% select(-t) %>% slice(1)


network.fort <- fortify(network)

# network.fort <- network.fort %>% select(-state) %>%
#   left_join(
#     Sim %>% select(-t) %>% slice(1) %>%
#       pivot_longer(cols = everything(), 
#                    names_to = "vertex.names", 
#                    values_to = "state") %>%
#       mutate(vertex.names = as.numeric(vertex.names),
#              state = factor(state, levels = c("S", "I", "R")))
#   )

# network.fort %>%
#   ggplot(aes(x, y, xend = xend, yend = yend)) +
#   geom_edges(color = "grey", alpha = 0.5) +
#   geom_nodes(aes(color = state), size = 2) + 
#   scale_color_manual(values = colors, label = c("S", "I", "R"), drop = FALSE) +
#   theme_blank()


Sim.long <- Sim %>%
  pivot_longer(cols = -t, 
               names_to = "vertex.names", 
               values_to = "state") %>%
  mutate(vertex.names = as.numeric(vertex.names),
         state = factor(state, levels = c("S", "I", "R"))) %>%
  mutate(time = format(t, digits = 3)) %>%
  # mutate(time = as.numeric(time))
  mutate(time = factor(str_remove(time, " "), ordered = T))


Sim.long.network <- Sim.long %>%
  left_join(network.fort)



p <- Sim.long.network %>%
  ggplot(aes(x, y, xend = xend, yend = yend)) +
  geom_edges(color = "grey", alpha = 0.5) +
  geom_nodes(aes(color = state), size = 2) + 
  scale_color_manual(values = colors, label = c("S", "I", "R"), drop = FALSE) +
  theme_blank()
p
# p <- p + transition_time(time) +
#   labs(title = "Time: {frame_time}")

# p <- p + transition_states(time, state_length = 0, wrap = FALSE) +
#   labs(title = "Time: {closest_state}")

p <- p +
  transition_manual(frames = factor(time, levels = unique(Sim.long.network$time))) +
  labs(title = "Time: {current_frame}")

nframes <- unique(Sim.long.network$time) %>% length()
nframes / 30

animate(p, nframes = nframes,
        renderer = gifski_renderer("../output/simulation1000_.gif"),
        fps = 25,
        end_pause = 8)

# animate(p, nframes = 2 * (unique(Sim.long.network$time) %>% length()),
#         renderer = gifski_renderer("../output/simulation.gif"))


# animate(p, renderer = ffmpeg_renderer("../output/simulation.mpg"))








# ---------------------------------------
#     Solving ODEs for SIR
# ---------------------------------------
N <- 5000
gamma <- 0.2
R0 <- 2.5

parameters <- c(beta = R0 * gamma,
                gamma = gamma, 
                N = N)
state <- c(S = N - 2, I = 2, R = 0)

SIR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- - beta / N * S * I
    dI <- beta / N * S * I - gamma * I
    dR <- gamma * I
    
    list(c(dS, dI, dR))
  })
}

t <- seq(0, 60, by = 0.01)

sim <- ode(y = state, times = t, func = SIR, parms = parameters) %>%
  as_tibble() %>% pivot_longer(-time, names_to = "compartment")

sim %>%
  ggplot() +
  geom_line(aes(x = time, y = value, group = compartment, color = factor(compartment, levels = c("S", "I", "R"))), size = 1) +
  scale_color_discrete(name = "State") +
  coord_cartesian(xlim = c(0, 60), y = c(0, N))



SimSummary.long %>%
  ggplot() +
  geom_line(aes(x = t, y = Count, group = interaction(iter, State), 
                color = factor(State, levels = c("S", "I", "R"))), alpha = 0.2) +
  geom_line(data = sim, 
            aes(x = time, y = value, group = compartment, color = factor(compartment, levels = c("S", "I", "R"))), size = 1) +
  scale_color_discrete(name = "State") +
  coord_cartesian(xlim = c(0, 60), y = c(0, N))



