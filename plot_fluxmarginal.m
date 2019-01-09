% plot marginals

function [maxp,exp_true,t,aux] = plot_fluxmarginal(x0,x1,mu,sigma,lb,ub,av,s,colormax)
    npoints = 4e5;
    exp_true = 0;
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    xlinsp = linspace(lb,ub, npoints);

    try
        t = truncate(pd, lb, ub);
        aux =  pdf(t, xlinsp);
        lbf = min(xlinsp(find(aux > 1e-5)));
        ubf = max(xlinsp(find(aux > 1e-5)));
    catch
        aux = xlinsp .* 0;
    if( sum(isnan(aux)) || sum(isinf(aux)) || nnz(aux) == 0 )
        xlinsp = linspace(x0,x1, npoints);
        pd = makedist('Exponential','mu',av);
        t = truncate(pd, x0, x1);
        exp_true = 1;
        aux =  pdf(t, xlinsp);
        lbf = min(xlinsp(find(aux > 1e-6)));
        ubf = max(xlinsp(find(aux > 1e-6)));
    end
    end
    plot(xlinsp, pdf(t, xlinsp),'Color',colormax);
    hold on
    area(xlinsp, pdf(t, xlinsp),'LineStyle','-','EdgeColor', colormax,'FaceColor',colormax,'FaceAlpha',0.3)
    maxp = max( pdf(t, xlinsp));
        
    
    
end
