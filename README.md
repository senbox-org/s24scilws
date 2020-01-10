S2-4Sci Land and Water Study Toolbox (s24scilws)
==========================

A toolbox for ...

Licenses
---------

...

Building from the source
------------------------------

Download and install the required build tools
	* Install J2SE 1.8 JDK and set JAVA_HOME accordingly. 
	* Install Maven and set MAVEN_HOME accordingly. 
	* Install git

Add $JAVA_HOME/bin, $MAVEN_HOME/bin to your PATH.

Clone the s24scilws source code and related repositories into a directory referred to as ${snap} from here on

    cd ${snap}
    git clone https://github.com/senbox-org/s24scilws.git
    git clone https://github.com/senbox-org/snap-desktop.git
    git clone https://github.com/senbox-org/snap-engine.git
    git clone https://github.com/senbox-org/s2tbx.git
    
Build SNAP-Engine:

    cd ${snap}/snap-engine
    mvn install

Build SNAP-Desktop:

    cd ${snap}/snap-desktop
    mvn install

Build s2tbx:

    cd ${snap}/s2tbx
    mvn install

Build s24scilws Toolbox:

    cd ${snap}/s24scilws
    mvn install
   
If unit tests are failing, you can use the following to skip the tests
   
    mvn clean
    mvn install -Dmaven.test.skip=true
	
Setting up IntelliJ IDEA
------------------------

1. Create an empty project with the ${snap} directory as project directory

2. Import the pom.xml files of snap-engine, snap-desktop and sen2cor as modules. Ensure **not** to enable
the option *Create module groups for multi-module Maven projects*. Everything can be default values.

3. Set the used SDK for the main project. A JDK 1.8 or later is needed.

4. Use the following configuration to run SNAP in the IDE:

    **Main class:** org.esa.snap.nbexec.Launcher
    **VM parameters:** -Dsun.awt.nopixfmt=true -Dsun.java2d.noddraw=true -Dsun.java2d.dpiaware=false
    All VM parameters are optional
    **Program arguments:**
    --userdir
    "${snap}/s24scilws/target/userdir"
    --clusters
    "${snap}/s2tbx/s2tbx-kit/target/netbeans_clusters/s2tbx;${snap}/s24scilws/s24scilws-kit/target/netbeans_clusters/s24scilws"
    --patches
    "${snap}/snap-engine/$/target/classes;${snap}/s2tbx/$/target/classes;${snap}/s24scilws/$/target/classes"
    **Working directory:** ${snap}/snap-desktop/snap-application/target/snap/
    **Use classpath of module:** snap-main

Enjoy!


