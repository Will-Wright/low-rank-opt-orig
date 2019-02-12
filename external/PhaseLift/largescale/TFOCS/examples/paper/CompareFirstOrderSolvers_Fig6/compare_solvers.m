%{
    Compares the various first-order algorithms applied to a smoothed Dantzig
    selector model with a known exact solution. Plots the relative error vs.
    applications of the linear operator A^*A.

    Should take a few minutes
%}

if exist('dtest_30dB_mu_25.mat','file');
    load dtest_30dB_mu_25.mat
    SS.type = '()'; SS.subs = { omega, ':' };
    N = length(xExact); M = length(b);
    Af = @(x) subsref(dct(x),SS);
    At = @(x) idct(subsasgn(zeros(N,size(x,2)),SS,x));
else
    disp('Regenerating the reference solution and problem data');
    opts          = [];
    opts.seed     = 83241;
    opts.smoothed = true;  % no continuation
    opts.SNR      = 30; 
    opts.type     = 'dct';
    opts.N        = 2^11;
    opts.M        = round(opts.N/4);
    [xExact,A,b,normA2,delta0,D,x0,xPlug,mu0,omega] = setup_Dantzig(opts);
    N = length(xExact); M = length(b);
    save('dtest_30dB_mu_25','D','Lexact','b','delta0','mu0','normA2','omega',...
        'xExact','xPlug');
end
% A = @(x,mode)linop_dantzig( M, N, Af, At, x, mode );
A = linop_handles( [M,N], Af, At );
x0 = xPlug;
z0 = zeros(N,1);

oExact   = norm(xExact,1) + mu0*norm(xExact-xPlug).^2/2;
nExact   = norm(xExact);
sExact   = xExact ~= 0;
err_norm = @(f,l,xk) norm( xk - xExact ) / nExact;
err_supp = @(f,l,xk) nnz( (xk~=0) ~= sExact );

opts        = [];
opts.tol    = 0;
opts.errFcn = { err_norm, err_supp };
opts.printEvery = 1000;
opts.maxcounts  = 10000;
opts.Lexact = Lexact;
solverList  = { 'AT', 'LLM', 'TS', 'N07', 'N83', 'GRA' };
legendStrings   = { 'fixed step','no restart','restart=100','restart=200','restart=400' };
betas       = [ 1, 0.5, 0.5, 0.5, 0.5 ]; % beta = 1 means no backtracking
restarts    = [ Inf, Inf, 100, 200, 400 ];
Lexacts     = [ Lexact, Inf, Inf, Inf, Inf ];
plot_data   = cell(3,length(solverList),length(betas));
tic
for si = 1 : length(solverList),
    opts.alg = solverList{si};
    
    % Only go through all the restart options for the AT solver
%     if isequal(opts.alg,'GRA'),
    if ~isequal(opts.alg,'AT'),
        nbetas = 2;
    else
        nbetas = length(betas);
    end
    for k = 1 : nbetas,
        opts.beta = betas(k);
        opts.restart = restarts(k);
        opts.Lexact = Lexacts(k);
        fprintf( '%s, beta=%g, restart=%d\n', solverList{si}, opts.beta, opts.restart );
        [ x2, out ] = solver_sDantzig( { A, D }, b, delta0, mu0, x0, z0, opts );
        plots{1,si,k} = out.counts(:,3);
        plots{2,si,k} = out.err(:,1);
        plots{3,si,k} = out.err(:,1) .* ~out.err(:,2);
    end
    
    if isequal(opts.alg,'AT')
        figure(si); clf;
        g = semilogy( plots{1:2,si,:}, 'linewidth', 1 );
        set(gca,'fontsize',16);
        legend(legendStrings{1:nbetas},'location','northeast');
        hold on
        for k = 1:nbetas,
            semilogy( plots{1,si,k}, plots{3,si,k}, 'linewidth', 3, 'color', get(g(k),'color') );
        end
        xlim([0,10000]);
        ylim([1e-9,1e-1]);
        xlabel('calls to $\mathcal{A}$ and $\mathcal{A}^*$','interpreter','latex','fontsize',16);
        ylabel('$\|x_k\mathrm{--}x_{\mu}^\star\|/\|x_\mu^\star\|$','interpreter','latex','fontsize',16);
        title(sprintf('%s method, various restart strategies',solverList{si}),'fontsize',12);
        %     print('-depsc2',sprintf('%s_restart.eps',solverList{si}));
    end
end
toc

%-- Make the figure that compares all solvers on one plot
figure(length(solverList)+1); clf;
g = semilogy( plots{1:2,:,2}, 'linewidth', 1 );
set(gca,'fontsize',12);
legend(solverList,'location','northeast');
hold on
for si = 1:length(solverList),
    semilogy( plots{1,si,2}, plots{3,si,2}, 'linewidth', 3, 'color', get(g(si),'color') );
    semilogy( plots{1,si,1}, plots{2,si,1}, 'linewidth', 1, 'linestyle', '--', 'color', get(g(si),'color') );
end
xlabel('calls to $\mathcal{A}$ and $\mathcal{A}^*$','interpreter','latex','fontsize',16);
ylabel('$\|x_k\mathrm{-}x_{\mu}^\star\|/\|x_\mu^\star\|$','interpreter','latex','fontsize',16);
title('All variants, fixed and backtracking steps, no restart','fontsize',20);
xlim([0,10000]);
ylim([1e-6,1e-1]);
% print -depsc2 combined.eps
%% print to PDF
figure(1);
lh = legend;
set(lh,'fontsize',18);
set(gca,'fontsize',16)
set(get(gca,'xlabel'),'fontsize',24)
set(get(gca,'ylabel'),'fontsize',24)
set(get(gca,'title'),'fontsize',20)

a = -.2;
orient landscape
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf,'PaperPosition',[a a (11-a) (8.5-a)]);
print('-dpdf','fig6b_recreation');

figure(length(solverList)+1)
lh = legend;
set(lh,'fontsize',18);
set(gca,'fontsize',16)
set(get(gca,'xlabel'),'fontsize',24)
set(get(gca,'ylabel'),'fontsize',24)
set(get(gca,'title'),'fontsize',20)

a = -.2;
orient landscape
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf,'PaperPosition',[a a (11-a) (8.5-a)]);
print('-dpdf','fig6a_recreation');



