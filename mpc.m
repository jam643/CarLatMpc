classdef mpc
    properties
        G
        H
        P_inv
        P
        f_1
    end
    methods
        function obj=mpc(A,B,C,N_pred,Q,R)
            [obj.H,obj.G] = obj.predictMats(A,B,C,N_pred);
            
            Q_rep = diag(repmat(diag(Q),N_pred,1));
            R_rep = diag(diag(R)*ones(N_pred,1));
            obj.P = obj.H'*Q_rep*obj.H+R_rep;
            obj.f_1 = obj.H'*Q_rep;
            obj.P_inv = inv(obj.P);
        end
        function [H,G] = predictMats(~,A,B,C,N_pred)
            q = size(C,1);
            n = size(A,2);

            G(1:q,1:n) = C*A;
            h(1:q,1) = C*B;
            H = zeros(N_pred*q,N_pred);
            for k = 2:N_pred
                h(end+1:end+q,1) = G(end-q+1:end,:)*B;
                G(end+1:end+q,:) = G(end-q+1:end,:)*A;
            end

            % lower triang matrix effect of U on prediction
            for k = 1:N_pred
                H(:,k) = [zeros(q*(k-1),1);h(1:end-(k-1)*q,1)];
            end
        end
        function U = solve(obj,x0,Yr)
            f = obj.f_1*(obj.G*x0 - Yr);
            U = -obj.P_inv*f;
        end
        function U = solveConstrained(obj,x0,Yr,A,B,Aeq,Beq,lb,ub)
            f = obj.f_1*(obj.G*x0 - Yr);
            U = quadprog(obj.P,f,A,B,Aeq,Beq,lb,ub);            
        end
    end
end
