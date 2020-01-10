package org.esa.s2tbx.s24scilws.topo;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.snap.core.util.ResourceInstaller;
import org.esa.snap.core.util.SystemUtils;
import org.esa.snap.runtime.Activator;

import java.io.IOException;
import java.nio.file.Path;

/**
 * Created by obarrile on 09/01/2020.
 */
public class TopoShadowActivator implements Activator {

    @Override
    public void start() {
        Path sourceDirPath = ResourceInstaller.findModuleCodeBasePath(getClass()).resolve("auxdata/s24scilws/topoShadow");
        Path auxdataDirectory = SystemUtils.getAuxDataPath().resolve("s2tbx/s24scilws/topoShadow");;
        if (auxdataDirectory == null) {
            SystemUtils.LOG.severe("TopoShadow configuration error: failed to retrieve auxdata path");
            return;
        }
        final ResourceInstaller resourceInstaller = new ResourceInstaller(sourceDirPath, auxdataDirectory);

        try {
            resourceInstaller.install(".*", ProgressMonitor.NULL);
        } catch (IOException e) {
            SystemUtils.LOG.severe("TopoShadow configuration error: failed to create " + auxdataDirectory);
            return;
        }
    }

    @Override
    public void stop() {
        // Purposely no-op
    }

    public static Path getAuxDataDir() {
        return SystemUtils.getAuxDataPath().resolve("s2tbx/s24scilws/topoShadow");
    }
}