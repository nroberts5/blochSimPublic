%% Initializing Parameter Values
params.beta = 1.0; 
params.M0 = 100.0;
params.T1W = 576e-3;
params.T1F = 280e-3;
params.R2s = 40;
params.psi = 25.0;
params.PDFF = 0.0;
params.rhoW = (1-params.PDFF)*params.M0;
params.rhoF = (params.PDFF)*params.M0;

%% Initializing Acquisition Parameters
acqParams.isPhantom = false;
acqParams.field_strength = 1.5;
acqParams.flip = 5.0;
acqParams.TE1 = 1.2e-3;
acqParams.DTE = 2.0e-3;
acqParams.necho = 6;
acqParams.tr = 15.2e-3;
acqParams.seed = 117.0;
acqParams.precessionDirection = 1.0;
acqParams.tes = linspace(acqParams.TE1,acqParams.TE1+(acqParams.necho-1)*acqParams.DTE,acqParams.necho);


%% Fat / Water Parameters
species(1).name = 'water'; % Water
species(1).frequency = [0];
species(1).relAmps = [1];
species(2).name = 'fat'; % Fat
species(2).frequency = -[3.80, 3.40, 2.60, 1.94, 0.39, -0.60];
species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
species(3).name = 'fatPhantom'; % Fat Phantom
species(3).frequency = -[3.80, 3.40, 2.60, 1.94, 0.39, -0.60]-0.1096;
species(3).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

if acqParams.isPhantom
    acqParams.ppm = species(3).frequency;
    acqParams.relAmps = species(3).relAmps;
else
    acqParams.ppm = species(2).frequency;
    acqParams.relAmps = species(2).relAmps;
end
acqParams.fp = 42.577 .* acqParams.ppm .* acqParams.field_strength;

%% Initializing Simulation Parameters
simParams.BlochSimReps = 1000; % # of repeated spins "per voxel" (needed to simulate spoiling)
simParams.blochNEX = 1000; % # number of TRs (i.e. number of times to repeat TR to achieve steady state)
simParams.forceSpoil = false;

%% Initialize the blochSim object
p = blochSim(); % init
p.addSpin(params.rhoW, params.T1W,1/params.R2s,params.psi,simParams.BlochSimReps); % add spin

% This loop is just for initializing Fat spins
for k=1:1:6
    m0 = params.rhoF*acqParams.relAmps(k);
    if m0>0
        p.addSpin(m0,params.T1F,1/params.R2s,params.psi+acqParams.fp(k),floor(simParams.BlochSimReps/6));
    end
end

%% Run the actual simulation
p = SGRE(p, params.beta, acqParams.flip, acqParams.tes, acqParams.necho, ...
    acqParams.tr, acqParams.seed, simParams.blochNEX, simParams.forceSpoil, 0);
signal_from_bloch = p.steadyStateSignal();

%% Prints results
fprintf('Steady State Magnitude:\n     %s\n',mat2str(round(abs(signal_from_bloch),4)));
fprintf('Steady State Phase:\n     %s\n',mat2str(round(angle(signal_from_bloch),4)));

%% Plotting with the Built in Tool
figure('Color','White','Name','Built in Plotting Tool'); clf;
subplot(1,2,1); p.plot_memory([],0,0); title('Magnitude');
subplot(1,2,2); p.plot_memory([],0,1); title('Phase');


%% Compare Result to Signal Equation
params.phi = 0;
signal_from_equation = SGRE_signal_eqution(params.beta, params.T1W, params.T1F, ...
    params.rhoW, params.rhoF, params.R2s, params.psi, params.phi,...
    acqParams.tes, deg2rad(acqParams.flip), acqParams.tr, acqParams.field_strength, acqParams.isPhantom);

figure('Color','White','Name','Comparing to Signal Model'); clf;
plotBloch(p, signal_from_equation);

%% Helper Functions

function p = SGRE(p, beta, tip, tes, nte, tr, PHI, NEX, forceSpoil, save0)
    k = 0;
    theta = 0;
    for nex=1:NEX
        % rf spoiling
        theta = mod(theta + k * PHI,360);
        % Increment k
        k = k+1;

        % Tip
        p.tip( beta*deg2rad(tip),deg2rad(theta));

        % Free Precess
        t = 0;
        for kt = 1:nte
            % Free Precess
            p.freeprecess(tes(kt)-t);
            % Measure
            p.save(save0+kt);
            t = tes(kt); 
        end

        % Decay
        p.freeprecess(tr-t);

        % Spoil
        if ~forceSpoil
            p.spoil();
        else
            p.forcespoil()
        end
    end
end

function out = SGRE_signal_eqution(beta, T1W, T1F, rhoW, rhoF, R2s, psi, phi, TEs, alpha, TR, field_strength, isPhantom)

species(1).name = 'water'; % Water
species(1).frequency = [0];
species(1).relAmps = [1]  ;
species(2).name = 'fat'; % Fat
species(2).frequency = -[3.80, 3.40, 2.60, 1.94, 0.39, -0.60];
species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
species(3).name = 'fatPhantom'; % Fat Phantom
species(3).frequency = -[3.80, 3.40, 2.60, 1.94, 0.39, -0.60]-0.1096;
species(3).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

if isPhantom
    ppm = species(3).frequency;
    relAmps = species(3).relAmps;
else
    ppm = species(2).frequency;
    relAmps = species(2).relAmps;
end
fp = 42.577 .* ppm .* field_strength;

T1_SEQ_TERM = @(alpha,TR,T1) (sin(beta.*alpha).*(1-exp(-TR./T1))) ./ (1-cos(beta.*alpha).*exp(-TR./T1));
T1W_TERM = T1_SEQ_TERM(alpha,TR,T1W);
T1F_TERM = T1_SEQ_TERM(alpha,TR,T1F);

out = zeros(size(TEs));
for te = TEs
    CN = sum(relAmps.*exp(1j.*2.*pi.*fp.*te));
    out(te==TEs) = (rhoW .* T1W_TERM + rhoF .* T1F_TERM .* CN) .* exp(-R2s .* te) .* exp(1j.* 2 .* pi .* te .* psi) .* exp(1j.*phi);
end

end

function plotBloch(p,equationSignal)

    subplot(1,2,1)
    p.plot_memory([],0,0);
    ax(1) = gca;
    lines1 = flip(ax(1).Children(arrayfun(@(x)strcmp(class(x),'matlab.graphics.chart.primitive.Line'),ax(1).Children)));
    for l = 1:numel(lines1)
        plot(get(ax(1),'XLim'),[abs(equationSignal(l)),abs(equationSignal(l))],'--','Color',lines1(l).Color)
    end

    subplot(1,2,2)
    p.plot_memory([],0,1);
    ax(2) = gca;
    lines1 = flip(ax(2).Children(arrayfun(@(x)strcmp(class(x),'matlab.graphics.chart.primitive.Line'),ax(2).Children)));
    for l = 1:numel(lines1)
        plot(get(ax(2),'XLim'),[angle(equationSignal(l)),angle(equationSignal(l))],'--','Color',lines1(l).Color)
    end

    linkaxes(ax,'x');
end