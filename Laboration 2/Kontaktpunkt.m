function contact_point = Kontaktpunkt(F, J, X)
    tol = 1e-10;
    iter = 0;
    maxiter = 100;
    normh = 1;
    hnew = 0;

    while normh > tol && iter < maxiter
        iter = iter + 1;
        h = -J(X(1), X(2), X(3))\F(X(1), X(2), X(3));
        X = X + h;
        normh = norm(h);
    end
    contact_point = X;
    % disp('iter/ x/ y/ a/ norm(h)');
    % disp([iter X(1) X(2) X(3) norm(h)])
end
