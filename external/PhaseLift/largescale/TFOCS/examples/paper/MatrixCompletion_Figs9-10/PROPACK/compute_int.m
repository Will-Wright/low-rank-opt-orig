function int = compute_int(mu,j,delta,eta,LL,strategy,extra)
%COMPUTE_INT:  Determine which Lanczos vectors to reorthogonalize against.
%
%      int = compute_int(mu,eta,LL,strategy,extra))
%
%   Strategy 0: Orthogonalize vectors v_{i-r-extra},...,v_{i},...v_{i+s+extra} 
%               with nu>eta, where v_{i} are the vectors with  mu>delta.
%   Strategy 1: Orthogonalize all vectors v_{r-extra},...,v_{s+extra} where
%               v_{r} is the first and v_{s} the last Lanczos vector with
%               mu > eta.
%   Strategy 2: Orthogonalize all vectors with mu > eta.
%
%   Notice: The first LL vectors are excluded since the new Lanczos
%   vector is already orthogonalized against them in the main iteration.

%  Rasmus Munk Larsen, DAIMI, 1998.

if (delta<eta)
    error('DELTA should satisfy DELTA >= ETA.')
end

if strategy == 0 && LL == 0 && extra == 0 && 3 == exist('computeInt','file')
    mu = abs(mu(1:j));
    
    I0 = find( mu >= delta);
    if isempty(I0)
        [mm_ignore,I0] = max(mu);
    end
%     int = zeros(j,1);
    
    % this is a mex file
    int = computeInt(I0,mu,eta);
%     for i = 1:length(I0)
        
        %             rList = I0(i):-1:1;
        %             r2 = find( abs(mu(rList))<eta | int(rList) == 1, 1 );
        
%         % Look backwards
%         for r=I0(i):-1:1
%             if mu(r)<eta || int(r)==1
%                 break;
%             else
%                 int(r) = 1;
%             end
%         end
%         
%         % Look forwards
%         for s=I0(i)+1:j
%             if mu(s)<eta || int(s)==1
%                 break;
%             else
%                 int(s) = 1;
%             end
%         end
    
%          int = computeInt(I0,mu,eta);


%     end
    int = find(int);
    int = int(:);
    return;
end

switch strategy
    case 0
        I0 = find(abs(mu(1:j))>=delta);
        if isempty(I0)
            [mm,I0] = max(abs(mu(1:j)));
        end
        int = zeros(j,1);
        for i = 1:length(I0)
            
            for r=I0(i):-1:1
                if abs(mu(r))<eta || int(r)==1
                    break;
                else
                    int(r) = 1;
                end
            end
            
            int(max(1,r-extra+1):r) = 1;
            for s=I0(i)+1:j
                if abs(mu(s))<eta || int(s)==1
                    break;
                else
                    int(s) = 1;
                end
            end
            int(s:min(j,s+extra-1)) = 1;
        end
        if LL>0
            int(1:LL) = 0;
        end
        int = find(int);
    case 1
        int=find(abs(mu(1:j))>eta);
        int = max(LL+1,min(int)-extra):min(max(int)+extra,j);
    case 2
        int=find(abs(mu(1:j))>=eta);
end
int = int(:);