function [particles, traj, landed_mask, veloc, meta] = move_particles(particles, dt, steps, g, recordTraj, sampleTimesSec)
% MOVE_PARTICLES  3D ballistic motion with sticky ground at y = 0.
% Semi-implicit Euler: v <- v + a*dt; x <- x + v*dt
% Ground plane at y=0: clamp y, zero velocity, stick thereafter.
%
% NEW: If you provide sampleTimesSec (seconds), only those frames are recorded.

    if nargin < 2 || isempty(dt),         dt = 0.01;   end
    if nargin < 3 || isempty(steps),      steps = 100; end
    if nargin < 4 || isempty(g),          g = 1.62;    end
    if nargin < 5 || isempty(recordTraj), recordTraj = false; end
    if nargin < 6,                         sampleTimesSec = []; end

    pos = particles.pos;  % (N x 3), [x y z]
    vel = particles.vel;  % (N x 3), [vx vy vz]
    N   = size(pos,1);

    % Ensure no particle starts below ground
    below = pos(:,2) < 0;
    if any(below)
        pos(below,2) = 0;
        vel(below,:) = 0;
    end
    active = true(N,1);    active(below) = false;

    a = [0, -g, 0];        % gravity in -Y

    % ------- build sampling indices (in frames) -------
    % Dense mode (legacy): record every step -> indices 1..steps+1
    denseMode = ~recordTraj || isempty(sampleTimesSec);

    if recordTraj
        if denseMode
            sample_idx = (1:(steps+1)).';                 % record every frame
            frame_times_sec = (0:steps).' * dt;
        else
            % Accept numeric vector or cell array of {struct(start,stop,step)}
            if iscell(sampleTimesSec)
                ts = [];
                for i = 1:numel(sampleTimesSec)
                    s = sampleTimesSec{i};
                    ts = [ts, s.start : s.step : s.stop]; %#ok<AGROW>
                end
            else
                ts = sampleTimesSec(:).';
            end
            % Keep only inside sim duration and include t=0 if missed
            t_end = steps * dt;
            ts = ts(ts >= 0 & ts <= t_end);
            if isempty(ts) || ts(1) > 0
                ts = [0, ts];            % ensure we capture the initial state
            end
            ts = unique(ts, 'stable');

            % Map seconds to frame indices (clamped)
            sample_idx = 1 + floor(ts./dt + 1e-9);
            sample_idx = min(steps+1, max(1, sample_idx(:)));
            % Remove duplicates after clamping
            [sample_idx, iu] = unique(sample_idx, 'stable');
            frame_times_sec = ((sample_idx-1) * dt).';     % actual recorded times
            ts = ts(iu);                                   %#ok<NASGU>  % kept for reference
        end

        % Preallocate recorded arrays to sampled length
        M = numel(sample_idx);
        traj = zeros(M, N, 3);
        landed_mask = false(M, N);
        veloc = zeros(M, N, 3);
    else
        traj = [];
        landed_mask = [];
        veloc = [];
        frame_times_sec = [];
        sample_idx = [];
    end

    % Write the very first recorded frame (t=0) if applicable
    if recordTraj
        write_ptr = 1;
        if sample_idx(write_ptr) == 1
            traj(write_ptr,:,:)   = pos;
            landed_mask(write_ptr,:) = ~active;
            veloc(write_ptr, :, :) = vel;
            write_ptr = write_ptr + 1;
        end
    end

    last_written_frame = 1;

    % -------- time integration ----------
    for s = 1:steps
        if any(active)
            % velocity update
            vel(active,1) = vel(active,1) + dt * a(1);
            vel(active,2) = vel(active,2) + dt * a(2);
            vel(active,3) = vel(active,3) + dt * a(3);
            % position update
            pos(active,1) = pos(active,1) + dt * vel(active,1);
            pos(active,2) = pos(active,2) + dt * vel(active,2);
            pos(active,3) = pos(active,3) + dt * vel(active,3);
        end

        % Ground collision
        hit = active & (pos(:,2) <= 0);
        if any(hit)
            pos(hit,2) = 0;
            vel(hit,:) = 0;
            active(hit) = false;
        end

        % Record if this frame is requested
        if recordTraj && write_ptr <= numel(sample_idx)
            current_frame = s + 1;  % frames are 1..steps+1
            while write_ptr <= numel(sample_idx) && sample_idx(write_ptr) == current_frame
                traj(write_ptr,:,:)      = pos;
                landed_mask(write_ptr,:) = ~active;
                veloc(write_ptr,:,:)     = vel;
                last_written_frame = current_frame;
                write_ptr = write_ptr + 1;
            end
        end

        % Early exit if all landed AND no later sample is requested
        if ~any(active) && (~recordTraj || write_ptr > numel(sample_idx))
            break
        end
    end

    % If we exited early but more samples were requested, pad with last state
    if recordTraj && write_ptr <= numel(sample_idx)
        for j = write_ptr:numel(sample_idx)
            traj(j,:,:)      = pos;
            landed_mask(j,:) = true;  % after everyone is stuck
            veloc(j,:,:)     = vel;
        end
    end

    % outputs
    particles.pos = pos;
    particles.vel = vel;

    if nargout >= 4
        meta = struct('frame_times_sec', frame_times_sec, ...
                      'dt', dt, 'steps', steps, ...
                      't_end_eff', (max(1, last_written_frame)-1)*dt);
    end
end