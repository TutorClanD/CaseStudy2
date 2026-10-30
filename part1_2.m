% ---------------------------------------------------------------------
function G = range_graph(D2, r)
% Adjacency matrix: link robots within distance r, no self-links.
G = double(D2 <= r^2);
G(logical(eye(size(G)))) = 0;
end

function tf = is_connected(G)
% True if every robot can reach every other robot through relay hops.
n = size(G, 1);
A = double((eye(n) + G) > 0);
R = A;
for k = 1:n
    if all(R(:))
        tf = true;
        return
    end
    Rnext = double((A*R) > 0);
    if isequal(Rnext, R)
        break
    end
    R = Rnext;
end
tf = all(R(:));
end
%  ---------------------------- Part I ---------------------------------
n = 100;
rng(105);
P0 = [5*rand(50,2); 5*rand(50,2) + [6 0]];

% Network building 
r = 2.1;
k = 11;
G1 = make_complete(n);
G2 = make_star(n);
G3 = make_ring(n,k); % when k=11, its edge count is closest to G4
G4 = make_disk(P0,r);

checkG4 = is_connected(G4);
% is_connected returns false when r=2.0; returns true when r=2.1

% Build Network Functions
function G = make_complete(n)
G = ones(n) - eye(n);
end

function G = make_star(n)
e1  = zeros(n,1);  
e1(1) = 1;   
one = ones(n,1);                

G = e1*one' + one*e1' - 2*(e1*e1');
end

function G = make_ring(n,k)
S = circshift(eye(n), 1);    
G = zeros(n);
Sj = eye(n);
for j = 1:k
    Sj = Sj*S;               
    G = G + Sj + Sj';
end
end

function G = make_disk(P,r)
n = size(P, 1);
s = zeros(n,1);
for i = 1:n
    s(i,1) = P(i,1)^2+P(i,2)^2;
end
delta = s*ones(1,n)+ones(n,1)*(s')-2*P*(P');
G = double(delta <= r^2);
G(logical(eye(n))) = 0;
end


% Analysis of network matrices
stats1 = graph_stats(G1);
stats2 = graph_stats(G2);
stats3 = graph_stats(G3);
stats4 = graph_stats(G4);

disp(stats3.edges);
disp(stats4.edges);




% function for calculating matrix info
function stats = graph_stats(G)
n = size(G, 1);

degree = G*ones(n,1);
meanDegree = mean(degree);
maxDegree = max(degree);
edges  = 0.5*ones(1,n)*G*ones(n,1);

A = eye(n) + G;
R = A;                      
diameter = Inf;            
for k = 1:n
    if all(R(:) > 0)       
        diameter = k;
        break
    end
    R = double(A*R > 0);     
end

stats = struct('degree', degree,'meanDegree',meanDegree, ...
    'maxDegree',maxDegree, 'edges', edges, 'diameter', diameter);
end


% Comparison Table
S = [stats1, stats2, stats3, stats4];      

Network    = ["T1"; "T2"; "T3"; "T4"];
Edges      = [S.edges]';
MeanDegree = [S.meanDegree]';
MaxDegree  = [S.maxDegree]';
Diameter   = [S.diameter]';

T = table(Network, Edges, MeanDegree, MaxDegree, Diameter)


% Graph plot
figure; 
hold on; 
axis equal; 
box on;

[X, Y] = gplot(G4, P0);            
plot(X, Y, '-', 'Color', [0.7 0.7 0.7]);                      

plot(P0(1:50,1),   P0(1:50,2),   'b.', 'MarkerSize', 14);    
plot(P0(51:100,1), P0(51:100,2), 'r.', 'MarkerSize', 14);     


xlabel('x'); ylabel('y');
title(sprintf('Disk Graph Network'));

%  ---------------------------- Part II ---------------------------------

% Setup: starting coordinates, networks, and the extra T4 with r = rmin + 0.4
x0 = P0(:,1);
y0 = P0(:,2);
r2 = round(r + 0.4, 1);                    % rmin + 0.4
G4b = make_disk(P0, r2);

nets    = {G1, G2, G3, G4, G4b};
netName = {'T1 complete', 'T2 star', 'T3 ring', 'T4 disk (rmin)', 'T4 disk (rmin+0.4)'};
nNet    = numel(nets);
K = 5000;                                  % long enough for every network to converge

% Sanity tests on small examples (5-robot star and ring)
for tcase = 1:2
    if tcase == 1
        Gt = make_star(5);  tname = '5-robot star';
    else
        Gt = make_ring(5,1); tname = '5-robot ring';
    end
    [WAt, WBt] = make_protocols(Gt);
    xt = [1; 2; 3; 4; 5];
    Xt = simulate_linear(WAt, xt, 10);
    fprintf('%s: size(X) = %s, final col == WA^10*x0 ? %d, WA*1==1 ? %d, WB*1==1 ? %d\n', ...
        tname, mat2str(size(Xt)), norm(Xt(:,end) - WAt^10*xt) < 1e-12, ...
        norm(WAt*ones(5,1) - 1) < 1e-12, norm(WBt*ones(5,1) - 1) < 1e-12);
end
% by eye: ring with k=1, A: x_1[1] = (x_5 + x_1 + x_2)/3 = (5+1+2)/3
[WAt, ~] = make_protocols(make_ring(5,1));
Xt = simulate_linear(WAt, [1;2;3;4;5], 1);
fprintf('ring check x_1[1] = %.4f (expected %.4f)\n', Xt(1,2), (5+1+2)/3);

% ---------------------------------------------------------------------
% - Prediction Log -
%
% E1 prediction: It is likely that the fastest acting set up will be the
%   Broadcast (Complete) simply due to the fact that every robot can hear
%   every other robot thus lowering travel delay. Contrasting this, T4 with
%   r_min will likely be the slowest since it connects to the fewest other
%   robots; therefore, it has the highest travel delay.
%
% E1 outcome   : We were right when predicting the fastest and the slowest
%   outcomes. There are no real surprises in the outcome because as the # involved 
%   with T increases the amount of robots connected to one another decreases, thus 
%   raising travel delay, which explains the decrease in overall speed. We
%   were also able to predict that r_min would be slower than r_min + 0.4,
%   since when we add the extra 0.4 to the minimum radius we increase the
%   amount of potential area that it can cover for communication thereby
%  decreasing the travel delay. So, the T4 with the smaller radius (r_min)
%  will be slower than the one with the larger radius (t_min + 0.4).
%
% E2 prediction: <predicted rendezvous point for each protocol; star distance from average>
%
% E2 outcome   : <were you right? explain any surprise>
% ---------------------------------------------------------------------

% ---------------------- E1: How fast? (Protocol A) -------------------
figure; hold on; box on;
kEps = zeros(nNet,1);
for t = 1:nNet
    [WA, ~] = make_protocols(nets{t});
    X = simulate_linear(WA, x0, K);
    Y = simulate_linear(WA, y0, K);
    dl = spread(X, Y);
    kEps(t) = rendezvous_time(dl);
    semilogy(0:K, dl, 'LineWidth', 1.5);
end
set(gca, 'YScale', 'log');
xlim([0 300]);                           
xlabel('step k'); ylabel('\delta[k]');
title('Protocol A: spread of the swarm'); legend(netName, 'Location', 'northeast');
yline(0.01*sqrt(sum((x0-mean(x0)).^2) + sum((y0-mean(y0)).^2)), 'k--', 'HandleVisibility', 'off');

Network = netName(:);
RendezvousTime_A = kEps;
E1 = table(Network, RendezvousTime_A)
[~, rankOrder] = sort(kEps);
disp('Fastest to slowest (Protocol A):'); disp(netName(rankOrder)');

% ---------------------- E2: Where do they meet? ----------------------
avgPt = [mean(x0), mean(y0)];
res = struct();
for t = 1:nNet
    Gt = nets{t};
    d  = Gt*ones(n,1);
    v  = d + 1;                            % v = (I+D)1
    [WA, WB] = make_protocols(Gt);

    % Predictions from the Checkpoint (c): conserved quantity / sum of weights
    predA = [v'*x0, v'*y0] / sum(v);       % degree-weighted average
    predB = avgPt;                         % plain average

    % Simulations
    XA = simulate_linear(WA, x0, K);  YA = simulate_linear(WA, y0, K);
    XB = simulate_linear(WB, x0, K);  YB = simulate_linear(WB, y0, K);
    simA = [mean(XA(:,end)), mean(YA(:,end))];
    simB = [mean(XB(:,end)), mean(YB(:,end))];

    res(t).name  = netName{t};
    res(t).errA  = norm(predA - simA);
    res(t).errB  = norm(predB - simB);
    res(t).distA = norm(predA - avgPt);    % how far Protocol A's meeting point is from the average
    res(t).distB = norm(predB - avgPt);    % always 0
    res(t).kA    = rendezvous_time(spread(XA, YA));
    res(t).kB    = rendezvous_time(spread(XB, YB));
    res(t).spreadA = spread(XA(:,end), YA(:,end));   % confirms all robots really met
    res(t).spreadB = spread(XB(:,end), YB(:,end));
end
Network      = netName(:);
PredErrA     = [res.errA]';
PredErrB     = [res.errB]';
DistFromAvgA = [res.distA]';
DistFromAvgB = [res.distB]';
TimeA        = [res.kA]';
TimeB        = [res.kB]';
E2 = table(Network, PredErrA, PredErrB, DistFromAvgA, DistFromAvgB, TimeA, TimeB)

% Star details
v = (G2*ones(n,1)) + 1;
fprintf('Star: Protocol A meets at (%.3f, %.3f); average is (%.3f, %.3f); distance = %.3f\n', ...
    (v'*x0)/sum(v), (v'*y0)/sum(v), avgPt(1), avgPt(2), res(2).distA);
fprintf('Star: rendezvous time A = %d, B = %d\n', res(2).kA, res(2).kB);

% ---------------------------------------------------------------------
% Local functions for Part II
% ---------------------------------------------------------------------
function [WA, WB] = make_protocols(G)
n = size(G,1);
d = G*ones(n,1);
D = diag(d);
I = eye(n);
L = D - G;
WA = (I + D) \ (I + G);                    % (I+D)^-1 (I+G)
WB = I - (1/max(d))*L;                     % eps = 1/max degree
end

function X = simulate_linear(W, x0, K)
n = numel(x0);
X = zeros(n, K+1);
X(:,1) = x0(:);
for k = 1:K
    X(:,k+1) = W*X(:,k);
end
end

function dl = spread(X, Y)
% delta[k] = sqrt(||C x[k]||^2 + ||C y[k]||^2), C = I - (1/n) 11'
n = size(X,1);
C = eye(n) - (1/n)*ones(n);
CX = C*X;
CY = C*Y;
dl = sqrt(sum(CX.^2, 1) + sum(CY.^2, 1));
end

function kE = rendezvous_time(dl)
% first step k with delta[k] <= 0.01*delta[0]  (dl(1) is delta[0])
idx = find(dl <= 0.01*dl(1), 1);
if isempty(idx)
    kE = Inf;
else
    kE = idx - 1;
end
end
