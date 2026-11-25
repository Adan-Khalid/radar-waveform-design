classdef TestRadarwave < matlab.unittest.TestCase
    %TESTRADARWAVE Unit tests of the radarwave package against closed-form theory.
    %   Run with: results = runtests('tests');  (from the matlab folder)

    methods (Test)
        function unitEnergy(tc)
            for s = {radarwave.lfm(64, 8), radarwave.phaseCode(radarwave.frankPhases(64)), ...
                    radarwave.ofdm(exp(1i * pi / 4 * (1:13)), 2)}
                tc.verifyEqual(sum(abs(s{1}).^2), 1, 'AbsTol', 1e-12);
            end
        end

        function barkerSidelobesExact(tc)
            for spc = [1 8]
                y = radarwave.correlate(radarwave.phaseCode(radarwave.barkerPhases(13), spc), ...
                    radarwave.phaseCode(radarwave.barkerPhases(13), spc));
                tc.verifyEqual(radarwave.pslDb(y), 20 * log10(1 / 13), 'AbsTol', 1e-9);
                tc.verifyEqual(radarwave.islDb(y), 10 * log10(12 / 169), 'AbsTol', 1e-9);
            end
        end

        function lfmApproachesSinc(tc)
            s = radarwave.lfm(1000, 8);
            y = radarwave.correlate(s, s);
            tc.verifyEqual(radarwave.pslDb(y), -13.26, 'AbsTol', 0.3);
            tc.verifyEqual(radarwave.islDb(y), -9.68, 'AbsTol', 0.1);
            tc.verifyEqual(radarwave.widthCells(y, -3, 8), 0.886, 'AbsTol', 0.01);
        end

        function echoPeakAtDelay(tc)
            s = radarwave.phaseCode(radarwave.p4Phases(32));
            x = zeros(200, 1);
            x(58:58+31) = s;
            [y, lag] = radarwave.correlate(x, s);
            [~, k] = max(abs(y));
            tc.verifyEqual(lag(k), 57);
        end

        function volumeIdentity(tc)
            s = radarwave.ofdm(exp(1i * pi / 2 * mod((0:12).^2, 4)), 2);
            K = 2 * numel(s);
            chi = radarwave.ambiguity(s, [], (0:K-1) * numel(s) / K);
            tc.verifyEqual(sum(abs(chi(:)).^2) / K, 1, 'RelTol', 1e-10);
        end

        function zeroDelayCutIsSampledSinc(tc)
            s = radarwave.lfm(32, 4);
            nu = linspace(-3, 3, 61);
            [chi, lags] = radarwave.ambiguity(s, [], nu);
            n = numel(s);
            expected = abs(sin(pi * nu) ./ (n * sin(pi * nu / n)));
            expected(nu == 0) = 1;
            tc.verifyEqual(abs(chi(lags == 0, :)), expected, 'AbsTol', 1e-12);
        end

        function lfmRidge(tc)
            s = radarwave.lfm(64, 8);
            for nut = [2 5 10]
                [chi, lags] = radarwave.ambiguity(s, [], nut);
                [pk, k] = max(abs(chi));
                tc.verifyEqual(pk, 1 - nut / 64, 'AbsTol', 0.01);
                tc.verifyEqual(lags(k) / 8, -nut, 'AbsTol', 1 / 8 + 1e-9);
            end
        end

        function taylorMatchesTaylorwin(tc)
            for n = [13 64 512]
                % compared after scaling to a peak of 1, so the test does not depend on
                % how a MATLAB release normalises taylorwin (losses and PSL do not either)
                w = radarwave.taylorWindow(n, 35, 5);
                ref = taylorwin(n, 5, -35);
                tc.verifyEqual(w / max(w), ref / max(ref), 'AbsTol', 1e-12);
                tc.verifyEqual(radarwave.windowFor('hann', n), hann(n), 'AbsTol', 1e-12);
                tc.verifyEqual(radarwave.windowFor('hamming', n), hamming(n), 'AbsTol', 1e-12);
            end
        end

        function windowLosses(tc)
            tc.verifyEqual(radarwave.windowLossDb(radarwave.windowFor('hamming', 4001)), 1.34, 'AbsTol', 0.01);
            tc.verifyEqual(radarwave.windowLossDb(radarwave.windowFor('hann', 4001)), 1.76, 'AbsTol', 0.01);
        end

        function quantizerLossAtMaxSteps(tc)
            expected = [1.9612 0.5495 0.1657 0.0504];
            for b = 1:4
                eta = radarwave.quantizerEfficiency(b, radarwave.gaussianOptimalFullScale(b));
                tc.verifyEqual(-10 * log10(eta), expected(b), 'AbsTol', 2e-3);
            end
            tc.verifyEqual(radarwave.quantizerEfficiency(1, 30), 2 / pi, 'RelTol', 1e-9);
        end

        function marcumMatchesNoncentralChiSquare(tc)
            snrDb = 11.749;
            pfa = 1e-4;
            viaNcx2 = 1 - ncx2cdf(-2 * log(pfa), 2, 2 * 10^(snrDb / 10));
            tc.verifyEqual(radarwave.pdMarcum(snrDb, pfa), viaNcx2, 'AbsTol', 1e-9);
            tc.verifyEqual(radarwave.requiredSnrDb(0.9, pfa), 11.749, 'AbsTol', 2e-3);
        end

        function lsFilterConstraintAndGain(tc)
            code = radarwave.phaseCode(radarwave.barkerPhases(13));
            ref = radarwave.lsMismatchedFilter(code, 39, 1e-9);
            y = radarwave.correlate(code, ref);
            tc.verifyEqual(max(abs(y)), 1, 'AbsTol', 1e-9);
            tc.verifyLessThan(radarwave.islDb(y), 10 * log10(12 / 169));
            loss = radarwave.mismatchLossDb(code, ref);
            tc.verifyGreaterThan(loss, 0);
            tc.verifyLessThan(loss, 1);
        end

        function midriseLevels(tc)
            for b = 1:4
                q = radarwave.quantizeReal(linspace(-3, 3, 100001), b, 1);
                lv = unique(q);
                tc.verifyEqual(numel(lv), 2^b);
                tc.verifyEqual(lv, -fliplr(lv), 'AbsTol', 1e-12);
            end
            q = radarwave.quantizeReal(linspace(-2, 2, 10001), 2, 1, 'midtread');
            tc.verifyEqual(unique(q), [-0.5 0 0.5], 'AbsTol', 1e-12);
        end

        function subtractiveDitherError(tc)
            rng(2, 'twister');
            x = (rand(400000, 1) - 0.5) * 1.6;
            u = rand(size(x)) - 0.5;
            e = radarwave.quantizeDither(x, 3, 1, 'midrise', 'subtractive', u) - x;
            d = 2 / 8;
            tc.verifyEqual(var(e), d^2 / 12, 'RelTol', 0.02);
            tc.verifyLessThan(abs(corr(e, x)), 0.01);
        end
    end
end
