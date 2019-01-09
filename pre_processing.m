function model = pre_processing(model0)

    model = [];
    Nfluxes = length(model0.lb);
    fprintf('Preprocessing..\n');
    ei = zeros(Nfluxes,1);
    lb = zeros(Nfluxes,1);
    ub = zeros(Nfluxes,1);
    options = optimset('linprog');
    options.Display = 'off';
    for i = 1:Nfluxes
        ei(i) = 1.0;
        sol = linprog(ei,[],[],model0.S,model0.b,model0.lb, model0.ub, [], options);
        lb(i) = max(model0.lb(i), sol(i));
        sol = linprog(-ei,[],[],model0.S,model0.b,model0.lb, model0.ub, [], options);
        ub(i) = min(model0.ub(i), sol(i));
        ei(i) = 0;
    end
    fixed = (lb == ub);
    fprintf('%d variables are fixed\n', nnz(fixed));
    notfixed = (lb < ub);
    model.S = model0.S(:,notfixed);
    model.b = model0.b - model0.S(:,fixed) * lb(fixed);
    model.lb = lb(notfixed);
    model.ub = ub(notfixed);
    model.rxns = model0.rxns(notfixed);
    model.rxnNames = model0.rxnNames(notfixed);
    model.mets = model0.mets;
    str = sprintf('%s pre-processed', model0.description);
    model.description = str; 


end