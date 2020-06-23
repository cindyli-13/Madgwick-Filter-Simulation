%% Madgwick Filter Implementation
%  Source: https://x-io.co.uk/open-source-imu-and-ahrs-algorithms/

classdef Madgwick < handle
    
    properties (Access = public)
        qEst = [1 0 0 0]    % orientation quaternion
        updateFreq = 512          % frequency at which filter is updated
        beta = 0.1                % represents mean zero gyro measurement error
        zeta = 0                  % gain used for gyro bias drift compensation
    end
    
    methods (Access = public)
        
        function obj = Madgwick(varargin)
            for i = 1:2:nargin
                if strcmp(varargin{i}, 'Quaternion'), obj.qEst = varargin{i+1};
                elseif strcmp(varargin{i}, 'UpdateFrequency'), obj.updateFreq = varargin{i+1};
                elseif strcmp(varargin{i}, 'Beta'), obj.beta = varargin{i+1};
                elseif strcmp(varargin{i}, 'Zeta'), obj.zeta = varargin{i+1};
                end
            end
        end
        
        % update() should be called at the filter's update frequency
        % gyro = [gx gy gz]
        % accel = [ax ay az]
        % mag = [mx my mz]
        function obj = update(obj, gyro, accel, mag)
            
            % put IMU measurements in 4D arrays
            gyro = [0 gyro];
            accel = [0 accel];
            mag = [0 mag];
            
            % for convenience
            q = obj.qEst;
            
            % angular rate quaternion
            qDot = 0.5 * [q(1)*gyro(1) - q(2)*gyro(2) - q(3)*gyro(3) - q(4)*gyro(4);
                          q(1)*gyro(2) + q(2)*gyro(1) + q(3)*gyro(4) - q(4)*gyro(3);
                          q(1)*gyro(3) - q(2)*gyro(4) + q(3)*gyro(1) + q(4)*gyro(2);
                          q(1)*gyro(4) + q(2)*gyro(3) - q(3)*gyro(2) + q(4)*gyro(1)]';
            
            % only compute feedback if accelerometer values are valid, (avoid Nan)
            if any(accel)
                % normalize accel measurement
                accel = accel / norm(accel);
                
                % objective function: fg(q, a) = q*g*q' - a (where g = [0 0 0 1])
                fg = 2 * [q(2)*q(4) - q(1)*q(3) - accel(2);
                          q(1)*q(2) + q(3)*q(4) - accel(3);
                          0.5 - q(2)^2 - q(3)^2 - accel(4)];
                
                % jacobian: Jg(q)
                jg = 2 * [-q(3),  q(2),   0;
                           q(4),  q(1),  -2*q(2);
                          -q(1),  q(4),  -2*q(3);
                           q(2),  q(3),   0];
                  
                % function gradient = Jg(q)fg(q, a)
                gradient = (jg * fg)';
                
                % only use magnetometer measurement if it is valid (avoid NaN)
                if any(mag)
                    % normalize mag measurement
                    mag = mag / norm(mag);
                    
                    % measured direction of Earth's magnetic field in earth frame
                    % h = q*mag*q' (quaternion rotation)
                    h = [0;
                         2*q(1)*(q(2)*mag(1) + q(3)*mag(4) - q(4)*mag(3)) + mag(2)*(q(1)^2 - q(2)^2 - q(3)^2 - q(4)^2);
                         2*(mag(2)*(q(1)*q(4) + q(2)*q(3)) + mag(4)*(q(3)*q(4) - q(1)*q(2))) + mag(3)*(q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2);
                         2*(mag(2)*(q(2)*q(4) - q(1)*q(3)) + mag(3)*(q(1)*q(2) + q(3)*q(4))) + mag(4)*(q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2)]';
                
                    % compensate for magnetic distortions
                    b = [0, sqrt(h(2)^2 + h(3)^2), 0, h(4)];
                    
                    % objective function: fb(q, b, m) = q*b*q' - m
                    fb = [2*b(2)*(0.5 - q(3)^2 - q(4)^2) + 2*b(4)*(q(2)*q(4) - q(1)*q(3)) - mag(2);
                          2*b(2)*(q(2)*q(3) - q(1)*q(4)) + 2*b(4)*(q(1)*q(2) + q(3)*q(4)) - mag(3);
                          2*b(2)*(q(1)*q(3) + q(2)*q(4)) + 2*b(4)*(0.5 - q(2)^2 - q(3)^2) - mag(4)];
                      
                    % jacobian: J(q, b)
                    jb = 2 * [-b(4)*q(3),                -b(2)*q(4) + b(4)*q(2),  b(2)*q(3);
                               b(4)*q(4),                 b(2)*q(3) + b(4)*q(1),  b(2)*q(4) - 2*b(4)*q(2);
                              -2*b(2)*q(3) - b(4)*q(1),   b(2)*q(2) + b(4)*q(4),  b(2)*q(1) - 2*b(4)*q(3);
                              -2*b(2)*q(4) + b(4)*q(2),  -b(2)*q(1) + b(4)*q(4),  b(2)*q(2)];
                
                    % function gradient = Jgb(q)fgb(q, a, b, m)
                    gradient = gradient + (jb * fb)';
                end
                
                % normalize gradient
                if any(gradient)
                    gradient = gradient / norm(gradient);
                end
                
                % compensate for gyro bias drift
                gyroErr = 2 * [q(1)*gradient(1) + q(2)*gradient(2) + q(3)*gradient(3) + q(4)*gradient(4);
                               q(1)*gradient(2) - q(2)*gradient(1) - q(3)*gradient(4) + q(4)*gradient(3);
                               q(1)*gradient(3) + q(2)*gradient(4) - q(3)*gradient(1) - q(4)*gradient(2);
                               q(1)*gradient(4) - q(2)*gradient(3) + q(3)*gradient(2) - q(4)*gradient(1)]';
                
                % integrate
                gyroBias = gyroErr / obj.updateFreq * obj.zeta;
                
                % compute compesnated gyro measurement
                gyro = gyro - gyroBias;
                
                % recompute angular rate quaternion
                qDot = 0.5 * [q(1)*gyro(1) - q(2)*gyro(2) - q(3)*gyro(3) - q(4)*gyro(4);
                              q(1)*gyro(2) + q(2)*gyro(1) + q(3)*gyro(4) - q(4)*gyro(3);
                              q(1)*gyro(3) - q(2)*gyro(4) + q(3)*gyro(1) + q(4)*gyro(2);
                              q(1)*gyro(4) + q(2)*gyro(3) - q(3)*gyro(2) + q(4)*gyro(1)]';
                
                % apply feedback
                qDot = qDot - obj.beta * gradient;
            end
            
            % integrate rate of change and add to current quaternion
            obj.qEst = obj.qEst + qDot / obj.updateFreq;
            if any(obj.qEst)
                obj.qEst = obj.qEst / norm(obj.qEst);
            end
        end
        
    end
    
end