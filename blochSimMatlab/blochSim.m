classdef blochSim < handle
    properties
            M0 
            T1 
            T2 
            df 
            repsPerPacket
            repsNormalizer
            M 
            
            spin_count
            PrecessionDirection
            
            time 
            tip_count 
            last_tip_time 

            mem 
            memkeys 

            exciterPhase 
            receiverPhase 
    end
    methods
        function obj = blochSim()
            obj.M0 = [];
            obj.T1 = [];
            obj.T2 = [];
            obj.df = [];
            obj.repsPerPacket = [];
            obj.repsNormalizer = [];
            obj.M = [];
            
            obj.spin_count = 0;
            obj.PrecessionDirection = 1;
            
            obj.time = 0;
            obj.tip_count = 0;
            obj.last_tip_time = 0;

            obj.mem = containers.Map('KeyType','single','ValueType','any');
            obj.memkeys = [];

            obj.exciterPhase = 0;
            obj.receiverPhase = 0;

    %         gyromagratio = 2*pi*42.58e6 %radians/T;
        end
        function newObj = deepCopy(obj)
            newObj = blochSim();
            
            newObj.M0 = obj.M0;
            newObj.T1 = obj.T1;
            newObj.T2 = obj.T2;
            newObj.df = obj.df;
            newObj.repsPerPacket = obj.repsPerPacket;
            newObj.repsNormalizer = obj.repsNormalizer;
            newObj.M = obj.M;
            
            newObj.spin_count = obj.spin_count;
            newObj.PrecessionDirection = obj.PrecessionDirection;
            
            newObj.time = obj.time;
            newObj.tip_count = obj.tip_count;
            newObj.last_tip_time = obj.last_tip_time;

            try
                newObj.mem = containers.Map(obj.mem.keys,obj.mem.values);
            end
            newObj.memkeys = obj.memkeys;

            newObj.exciterPhase = obj.exciterPhase;
            newObj.receiverPhase = obj.receiverPhase;
        end
        function reset(obj)
            obj.freeprecess(1e25)
            obj.time = 0;
            obj.tip_count = 0;
            obj.last_tip_time = 0;
        end
        function addSpin(obj,M0,T1,T2,df,reps)
            obj.M0 = [obj.M0, M0];
            obj.T1 = [obj.T1, T1];
            obj.T2 = [obj.T2, T2];
            obj.df = [obj.df, df];
            obj.repsPerPacket = [obj.repsPerPacket, reps];
            
            for rep = 1:reps
                obj.M = horzcat(obj.M,[0;0;M0]);
                obj.repsNormalizer = horzcat(obj.repsNormalizer,reps);
            end
            obj.spin_count = obj.spin_count + 1;
        end
        function Rx = xrot(obj,theta)
            Rx = [1,0,0;0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];
                end
        function Ry = yrot(obj,theta)
            Ry = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
        end
        function Rz = zrot(obj,theta)
            Rz = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
        end
        function Rth=throt(obj,phi,theta)
              Rz = obj.zrot(-theta);
              Ry = obj.yrot(phi);
              Rth = (Rz\Ry)*Rz;
        end 
        function tip(obj, alpha, phase)
            obj.exciterPhase = phase;
            obj.receiverPhase = phase;
            
            obj.M = obj.throt(alpha,obj.exciterPhase)*obj.M;
            obj.tip_count = obj.tip_count+1;
            
            obj.last_tip_time = obj.time;
        end
        function setReceiverPhase(obj, phase)
            obj.receiverPhase = phase;
        end
        function spoil(obj,low_ext,up_ext)
            if nargin < 2
                low_ext = -2*pi;
            end
            if nargin < 3
                up_ext = 2*pi;
            end
            %% Spoiling Option 2 (Each species recieves the same spoiling)
            for packet = 1:numel(obj.repsPerPacket)
                M_packet = obj.M(:, sum(obj.repsPerPacket(1:packet-1))+1 : sum(obj.repsPerPacket(1:packet)));
                phases = linspace(low_ext, up_ext, size(M_packet,2));
                for k=1:size(M_packet,2)
                    M_packet(:,k) = obj.zrot(phases(k))*M_packet(:,k);
                end
                obj.M(:, sum(obj.repsPerPacket(1:packet-1))+1 : sum(obj.repsPerPacket(1:packet))) = M_packet;
            end
        end
        function forcespoil(obj)
            obj.M(1:2,:) = 0;
        end
        function [A,B] = AB(obj, T, packet)
            phi = obj.PrecessionDirection*2*pi*obj.df(packet) * T;
            E1 = exp(-T/obj.T1(packet));
            E2 = exp(-T/obj.T2(packet));
            
            A = [E2,0,0;0,E2,0;0,0,E1] * obj.zrot(phi);
            B = [0;0;obj.M0(packet)*(1-E1)];
        end
        function freeprecess(obj,time_step)
            for packet = 1:numel(obj.repsPerPacket)
                M_packet = obj.M(:, sum(obj.repsPerPacket(1:packet-1))+1 : sum(obj.repsPerPacket(1:packet)));
                [Afp, Bfp] = obj.AB(time_step,packet);
                M_packet = Afp*M_packet + Bfp;
                obj.M(:, sum(obj.repsPerPacket(1:packet-1))+1 : sum(obj.repsPerPacket(1:packet))) = M_packet;
            end
            obj.time = obj.time + time_step;
        end
        function sig = transverseSignal(obj)
            sig = sum((obj.M(1,:) + 1j*obj.M(2,:)) * exp(-1j*obj.receiverPhase) ./ obj.repsNormalizer,2);
        end
        function sig = longitudinalSignal(obj)
            sig = sum(obj.M(3,:) ./ obj.repsNormalizer,2);
        end
        function save(obj, mem_key)
            if ~ismember(mem_key,obj.memkeys)
                obj.mem(mem_key) = [];
                obj.memkeys = horzcat(obj.memkeys,mem_key);
            end
            obj.mem(mem_key) = vertcat(obj.mem(mem_key), [obj.transverseSignal(), obj.longitudinalSignal(), obj.time]);
        end
        function [trans_signal, long_signal, time_pts] = download_memory(obj,mem_key)
               data = obj.mem(mem_key);
               trans_signal = data(:,1).';
               long_signal  = data(:,2).';
               time_pts     = abs(data(:,3).');
        end
        function sig_ss = steadyStateSignal(obj, memkeys_to_record, lastNpts)
            if nargin<2
                memkeys_to_record = obj.memkeys;
            elseif isempty(memkeys_to_record)
                memkeys_to_record = obj.memkeys;
            end
            
            if nargin<3
                lastNpts = 16;
            end
            
            sig_ss = zeros(1,numel(memkeys_to_record));
            for memkey = memkeys_to_record
                sig = obj.download_memory(memkey);
                strt = max([numel(sig)-lastNpts+1,1]);
                sig_ss(memkey == memkeys_to_record) = median(sig(strt:end));
            end
        end
        function plot_handles = plot_memory(obj, memkeys_to_plot, plot_long, plot_angle)
            if nargin<2
                memkeys_to_plot = obj.memkeys;
            elseif isempty(memkeys_to_plot)
                memkeys_to_plot = obj.memkeys;
            end
            if nargin<3
                plot_long = 0;
            end
            if nargin<4
                plot_angle = 0;
            end
            place_text = @(x, yfunc, prct, txt, color) ...
                text(prctile(x,prct),yfunc(prctile(x,prct)), txt, 'Color', color, 'FontSize', 18);
            plot_handles = [];
            for memkey = memkeys_to_plot
               [trans_signal, long_signal, time_pts] = obj.download_memory(memkey);
               if plot_angle
                   yfuncTRANS = griddedInterpolant(time_pts,angle(trans_signal));
               else
                   yfuncTRANS = griddedInterpolant(time_pts,abs(trans_signal));
               end
               trans_plt = plot(time_pts,yfuncTRANS(time_pts),'-');  
               place_text(time_pts, yfuncTRANS, mod(20+10*(memkey-1),90), sprintf('S_{%g}^{T}',memkey), trans_plt.Color)
               hold on;
               if plot_long
                  if plot_angle
                       yfuncLONG = griddedInterpolant(time_pts,angle(long_signal)); 
                   else
                       yfuncLONG = griddedInterpolant(time_pts,abs(long_signal)); 
                   end          
                   long_plt = plot(time_pts,yfuncLONG(time_pts),'--','Color',trans_plt.Color);
                   place_text(time_pts, yfuncLONG , mod(20+10*(memkey-1),90), sprintf('S_{%g}^{L}',memkey), trans_plt.Color)
               end
               plot_handles = [plot_handles,trans_plt];
            end
            
            xlabel('Time (s)')
            if plot_angle
                ylabel('Signal Phase')
            else
                ylabel('Signal Magnitude')
            end
        end
        function Mplot(obj,plot_all_M,color1,color2)
           if nargin < 2
                plot_all_M = 0;
            end
            if nargin < 3
                color1 = 'k';
            end
           if nargin < 4
                color2 = color1;
            end
            if plot_all_M
                for m = 1:size(obj.M,2)
                    hold on;
                    vectarrow([0,0,0],obj.M(:,m).',color1,1)
                end
            end
            hold on;
            vectarrow([0,0,0],mean(obj.M,2).',color2,5)
        end
    end
end
