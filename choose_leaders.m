function S = choose_leaders(G, P0, faulty, g, m)
%CHOOSE_LEADERS  Choose which robots become leaders (ESE 105, Case Study II, Part III).
%
%   S = choose_leaders(G, P0, faulty, g, m)
%
%   Inputs
%     G       n-by-n adjacency matrix of the communication network
%     P0      n-by-2 landing positions; row i is robot i's (x, y)
%     faulty  vector of indices of robots with broken drives
%     g       1-by-2 rally point
%     m       maximum number of leaders
%
%   Output
%     S       vector of at most m distinct robot indices, none of them faulty
%
%   Rules (see Section D.4 of the handout)
%     - Use only the network's structure: anything you can compute from G,
%       P0, faulty and g with the tools from Part I (degrees, hop counts,
%       distances, products and powers of G, ...).
%     - Do NOT simulate the swarm, and do NOT call score_mission.
%     - Must run in under 2 seconds on one deployment.
%     - Do not change the function name or the order of the inputs.
%       Helper functions are allowed, but put them in this file, below
%       this function, so that this single file is all we need to run it.
%
%   Team: <CS2 Group 9>
%   Members: <Declan Fowler, Jasper James, & Jimmy, Xu>
%
%   Strategy (2-3 sentences describing your rule, in words):
%   Each robot j is given a "need" equal to its squared distance from g (far
%   robots have the most to gain), discounted if its neighborhood is mostly
%   faulty robots. A candidate leader i is valued by how much of that need
%   it can reach: the sum over j of need_j times K_ji, where K = average of
%   W^1..W^H is the few-hop averaging kernel built from G (W = (I+D)^-1 (I+G)).
%   We pick leaders greedily, and after each pick we reduce the remaining
%   need of every robot that leader already covers, so the m leaders end up
%   spread over the different clusters instead of piled in one place.

n = size(G, 1);
working = true(n, 1);
working(faulty) = false;
m = min(m, nnz(working));

H = 4;                                   % number of hops of influence considered (Needs to be Tested and changed)

% ---- structure from G (Part I tools: degrees, products/powers of G) ----
d = G*ones(n, 1);                        % degree vector
W = (eye(n) + G) ./ (d + 1);             % (I + D)^-1 (I + G): row i scaled by 1/(1+d_i)
K = zeros(n);
Wh = eye(n);
for h = 1:H
    Wh = Wh*W;                           % W^h: where an h-hop average of j's neighborhood comes from
    K = K + Wh/H;                        % K(j,i): how strongly robot i influences robot j
end

% ---- need: distance^2 to g, discounted when j mostly listens to faulty robots ----
dist2 = sum((P0 - g).^2, 2);             % squared distance of each robot from g
toFaulty = min(1, sum(K(:, faulty), 2)); % share of j's influence that comes from faulty robots
need = dist2 .* (1 - toFaulty) .* working;

% ---- greedy selection with diminishing returns ----
Kn = K / max(K(:));                      % influence scaled into [0,1] so coverage is a fraction
resid = ones(n, 1);                      % fraction of each robot's need not yet covered
S = zeros(1, m);
for t = 1:m
    gain = ((need .* resid)' * K)';      % value of making each robot a leader
    gain(~working) = -Inf;
    gain(S(1:t-1)) = -Inf;
    [~, i] = max(gain);
    S(t) = i;
    resid = resid .* (1 - min(1, Kn(:, i)));
end

end
