%% Case Study II, Part III: Mission Design and Competition
% Uses make_deployment.m and choose_leaders.m from this folder.
clear; clc; close all;

lambdaLeader = 0.5;   % leader responsiveness
T = 150;              % score horizon
m = 4;                % leader budget

%% D.1  Test simulate_anchored on a tiny example (5-robot ring)
% Hand check: with all lambda = 1 it must reduce to plain W*x.
n5 = 5;  G5 = double(abs((1:n5)' - (1:n5)) == 1 | abs((1:n5)' - (1:n5)) == n5-1);
W5 = (eye(n5) + G5) ./ (G5*ones(n5,1) + 1);
x5 = [1;2;3;4;5];
X5 = simulate_anchored(W5, ones(n5,1), zeros(n5,1), x5, 3);
fprintf('lambda=1 matches W^3*x0 ? %d\n', norm(X5(:,end) - W5^3*x5) < 1e-12);
% lambda = 0 for everyone: nobody moves, state stays at the anchor.
X5 = simulate_anchored(W5, zeros(n5,1), x5, x5, 3);
fprintf('lambda=0 stays put ? %d\n', norm(X5(:,end) - x5) < 1e-12);

%% D.1  Steady state on a practice deployment, no leaders
D = make_deployment(1);
fprintf('make_disk reproduces D.G ? %d\n', isequal(make_disk(D.P0, D.r), D.G));
[lam, a, b] = build_setup(D, []);               % no leaders
W = protocol_A(D.G);
n = size(D.G, 1);
Lam = diag(lam);
% x* = Lam W x* + (I-Lam) a   =>   (I - Lam W) x* = (I - Lam) a
xs = (eye(n) - Lam*W) \ ((eye(n) - Lam)*a);
ys = (eye(n) - Lam*W) \ ((eye(n) - Lam)*b);
X = simulate_anchored(W, lam, a, D.P0(:,1), 5000);
Y = simulate_anchored(W, lam, b, D.P0(:,2), 5000);
fprintf('steady state vs long run: %.2e\n', norm([xs ys] - [X(:,end) Y(:,end)], 'fro'));
fprintf('no leaders: mean steady-state position = (%.2f, %.2f);  mean faulty position = (%.2f, %.2f)\n', ...
    mean(xs), mean(ys), mean(D.P0(D.faulty,1)), mean(D.P0(D.faulty,2)));
% Why: faulty robots never move, so every steady state is a weighted mix of
% their landing positions; the whole swarm drifts to the faulty zone.

%% D.5 Step 1  Three simple strategies vs. mine on practice deployments
seeds = 1:20;
names = {'Random', 'Highest degree', 'Closest to g', 'Mine'};
J = zeros(numel(seeds), 4);
for q = 1:numel(seeds)
    Dq = make_deployment(seeds(q));
    rng(seeds(q));
    S = {rule_random(Dq, m), rule_degree(Dq, m), rule_closest(Dq, m), ...
         choose_leaders(Dq.G, Dq.P0, Dq.faulty, Dq.g, m)};
    for s = 1:4
        J(q, s) = score_mission(Dq, S{s}, lambdaLeader, T, m);
    end
end

fprintf('\n%-16s %8s %8s %8s\n', 'Strategy', 'Mean', 'Worst', 'Best');
for s = 1:4
    fprintf('%-16s %8.2f %8.2f %8.2f\n', names{s}, mean(J(:,s)), min(J(:,s)), max(J(:,s)));
end

figure; boxplot(J, 'Labels', names); grid on;
ylabel('J (% reduction in distance to g)');
title(sprintf('Score distributions over %d practice deployments', numel(seeds)));

% Timing check for D.4 (must be < 2 s)
Dt = make_deployment(7);
tic; choose_leaders(Dt.G, Dt.P0, Dt.faulty, Dt.g, m); fprintf('choose_leaders time: %.3f s\n', toc);

%% Look closely at one failure: worst deployment for "closest to g"
[~, qf] = min(J(:,3));
Df = make_deployment(seeds(qf));
Sc = rule_closest(Df, m);  Sm = choose_leaders(Df.G, Df.P0, Df.faulty, Df.g, m);
figure; hold on; axis equal; box on;
[gx, gy] = gplot(Df.G, Df.P0);  plot(gx, gy, '-', 'Color', [0.85 0.85 0.85]);
plot(Df.P0(:,1), Df.P0(:,2), 'k.', 'MarkerSize', 10);
plot(Df.P0(Df.faulty,1), Df.P0(Df.faulty,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(Df.P0(Sc,1), Df.P0(Sc,2), 'bs', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(Df.P0(Sm,1), Df.P0(Sm,2), 'go', 'MarkerSize', 12, 'LineWidth', 2);
plot(Df.g(1), Df.g(2), 'm*', 'MarkerSize', 14, 'LineWidth', 2);
legend({'links','robots','faulty','closest-to-g leaders','my leaders','rally point g'});
title(sprintf('Seed %d: J(closest) = %.1f, J(mine) = %.1f', seeds(qf), J(qf,3), J(qf,4)));

%% ---------------------------- Local functions ----------------------------
function W = protocol_A(G)
n = size(G, 1);
W = (eye(n) + G) ./ (G*ones(n,1) + 1);          % (I + D)^-1 (I + G)
end

function X = simulate_anchored(W, lambda, a, x0, K)
% Column k+1 of X is x[k] for x[k+1] = Lam*W*x[k] + (I - Lam)*a.
n = numel(x0);
X = zeros(n, K+1);
X(:,1) = x0(:);
lambda = lambda(:);  a = a(:);
for k = 1:K
    X(:,k+1) = lambda .* (W*X(:,k)) + (1 - lambda) .* a;
end
end

function [lam, a, b] = build_setup(D, S, lambdaLeader)
% lambda and anchors for faulty robots and leaders S (D.1).
if nargin < 3, lambdaLeader = 0.5; end
n = size(D.G, 1);
lam = ones(n, 1);   a = D.P0(:,1);   b = D.P0(:,2);
lam(D.faulty) = 0;                               % faulty: anchor = landing spot
lam(S) = lambdaLeader;                           % leaders: anchor = g
a(S) = D.g(1);   b(S) = D.g(2);
end

function J = score_mission(D, S, lambdaLeader, T, m)
% Average percent reduction in mean squared distance to g (D.3).
S = S(:);
n = size(D.G, 1);
working = true(n, 1);  working(D.faulty) = false;
assert(numel(S) <= m && numel(unique(S)) == numel(S), 'S must have at most m distinct robots.');
assert(all(working(S)), 'S contains a faulty robot.');
assert(all(S >= 1 & S <= n & S == round(S)), 'S has invalid indices.');

W = protocol_A(D.G);
e = zeros(T, 2);                                 % columns: baseline, new
for c = 1:2
    if c == 1, Sc = []; else, Sc = S; end
    [lam, a, b] = build_setup(D, Sc, lambdaLeader);
    X = simulate_anchored(W, lam, a, D.P0(:,1), T);
    Y = simulate_anchored(W, lam, b, D.P0(:,2), T);
    X = X(:, 2:end);  Y = Y(:, 2:end);           % steps k = 1..T
    dist2 = (X(working,:) - D.g(1)).^2 + (Y(working,:) - D.g(2)).^2;
    e(:,c) = mean(dist2, 1)';
end
J = (100/T) * sum(1 - e(:,2) ./ e(:,1));
end

function S = rule_random(D, m)
w = setdiff(1:size(D.G,1), D.faulty);
S = w(randperm(numel(w), m));
end

function S = rule_degree(D, m)
d = D.G * ones(size(D.G,1), 1);
d(D.faulty) = -Inf;
[~, idx] = sort(d, 'descend');
S = idx(1:m)';
end

function S = rule_closest(D, m)
dist2 = sum((D.P0 - D.g).^2, 2);
dist2(D.faulty) = Inf;
[~, idx] = sort(dist2, 'ascend');
S = idx(1:m)';
end
