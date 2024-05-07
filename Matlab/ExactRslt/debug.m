clear; clc;

alternative_psi = rand(1, 4);
config = [1, 1, 0, 1];
site1 = 1; site2 = 2;

if config(site1) == config(site2)
    for config1 = 0:1
        for config2 = 0:1
            if config1 == config2
                config = config1 * 2 + config2
            end
        end
    end
else
    for config1 = 0:1
        for config2 = 0:1
            if config1 ~= config2
                config = config1 * 2 + config2
            end
        end
    end
end