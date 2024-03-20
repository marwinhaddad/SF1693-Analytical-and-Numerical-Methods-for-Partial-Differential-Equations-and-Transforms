clc;
format 

% Test problem 1
% rng(610729)
% Test problem 2
% rng(YYMMDD) Replace YYMMDD by first group member’s date of birth
% rng(980520)
% Test problem 3 applicable only if two persons in the group
% Test problem 3
% rng(YYMMDD) Replace YYMMDD by second member’s date of birth



% m=5;
% n=12;
% A=[randi([0 m],m,n-m) eye(m)];
% b=randi([m 2*m],m,1);
% c=[-randi([1 n-m],n-m,1) ; zeros(m,1) ];
% Basis = [n-m+1:n];

A = [3 2 1 1 0 0; % protein
     4 2 1 0 1 0; % grönsaker
     3 1 3 0 0 1]; % ris

b = [500; 600; 800];

[m, n] = size(A);

c = [-70; -40; -30; zeros(m, 1)];

% A = [3 2 1 3 3 1;
%      2 4 2 1 2 1;
%      1 2 3 2 3 3];
%  
%  b = [14; 16; 10];
%  
%  c = [2 3 2 2 3 2]';
%  
%  [m, n] = size(A);
 
Basis = [n-m+1:n];



xb = simplex(A, b, c, Basis)

function xb = simplex(A, b, c, Basis)
    % create vBasis
    vBasis = [];
    for n = 1:size(A, 2)
        if ~ismember(n, Basis)
            vBasis(end+1) = n;
        end
    end
    
    % initalize solution vector and iterations
    iter = 0;
    xb = zeros(size(A, 2), 1);
    while iter < 10
        % determine Ab, Av, cb and cv...
        Ab = A(:, Basis);
        Av = A(:, vBasis);

        cb = c(Basis);
        cv = c(vBasis);

        % ... to calculate bbar, y, rv
        bBar = Ab\b
        y = Ab'\cb;
        
        rv = cv - Av' * y
        if all(rv >= 0)
            disp("Optimal solution found")
            xb(Basis) = bBar;
            return
        end

        % choose q and calculate a_vqBar
        [~, q] = min(rv);
        vq = vBasis(q);

        a_vq = A(:, vq);
        a_vqBar = Ab\a_vq
        
        if all(a_vqBar <= 0)
            disp("solution does not exist")
            return
        end

        % choose p and interchange vBasis(q) and Basis(p)
        tmax = inf;
        p = 0;
        for i = 1:length(bBar)
            if a_vqBar(i) > 0
                ttemp = bBar(i) / a_vqBar(i);
                if ttemp < tmax
                    tmax = ttemp;
                    p = i;
                end
            end
        end
        
        vBasis(q) = Basis(p);
        Basis(p) = vq;
        
        iter = iter + 1
    end
    disp("Maximum amount of iterations exceeded")
    xb = NaN;
end
