// $Id: Util.java,v 1.29 2009/11/24 23:56:35 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URI;
import java.net.URL;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Locale;
import java.util.Set;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.JApplet;
import javax.swing.JOptionPane;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.tol.TolParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterUtil;

final class Util {

    private final static String[] AVAILABLE_FONT_FAMILIES_SORTED = GraphicsEnvironment.getLocalGraphicsEnvironment()
                                                                         .getAvailableFontFamilyNames();
    static {
        Arrays.sort( AVAILABLE_FONT_FAMILIES_SORTED );
    }

    final static void addPhylogeniesToTabs( final Phylogeny[] phys,
                                            final String default_name,
                                            final String full_path,
                                            final Configuration configuration,
                                            final MainPanel main_panel ) {
        if ( phys.length > Constants.MAX_TREES_TO_LOAD ) {
            JOptionPane.showMessageDialog( main_panel, "Attempt to load " + phys.length
                    + " phylogenies,\ngoing to load only the first " + Constants.MAX_TREES_TO_LOAD, Constants.PRG_NAME
                    + " more than " + Constants.MAX_TREES_TO_LOAD + " phylogenies", JOptionPane.WARNING_MESSAGE );
        }
        int i = 1;
        for( final Phylogeny phy : phys ) {
            if ( i <= Constants.MAX_TREES_TO_LOAD ) {
                String my_name = "";
                String my_name_for_file = "";
                if ( phys.length > 1 ) {
                    if ( !ForesterUtil.isEmpty( default_name ) ) {
                        my_name = new String( default_name );
                    }
                    if ( !ForesterUtil.isEmpty( full_path ) ) {
                        my_name_for_file = new String( full_path );
                    }
                    else if ( !ForesterUtil.isEmpty( default_name ) ) {
                        my_name_for_file = new String( default_name );
                    }
                    String suffix = "";
                    if ( my_name_for_file.indexOf( '.' ) > 0 ) {
                        suffix = my_name_for_file.substring( my_name_for_file.lastIndexOf( '.' ), my_name_for_file
                                .length() );
                        my_name_for_file = my_name_for_file.substring( 0, my_name_for_file.lastIndexOf( '.' ) );
                    }
                    if ( !ForesterUtil.isEmpty( my_name_for_file ) ) {
                        my_name_for_file += "_";
                    }
                    if ( !ForesterUtil.isEmpty( phy.getName() ) ) {
                        my_name_for_file += phy.getName().replaceAll( " ", "_" );
                    }
                    else if ( phy.getIdentifier() != null ) {
                        final StringBuffer sb = new StringBuffer();
                        if ( !ForesterUtil.isEmpty( phy.getIdentifier().getProvider() ) ) {
                            sb.append( phy.getIdentifier().getProvider() );
                            sb.append( "_" );
                        }
                        sb.append( phy.getIdentifier().getValue() );
                        my_name_for_file += sb;
                    }
                    else {
                        my_name_for_file += i;
                    }
                    if ( !ForesterUtil.isEmpty( my_name ) && ForesterUtil.isEmpty( phy.getName() )
                            && ( phy.getIdentifier() == null ) ) {
                        my_name = my_name + " [" + i + "]";
                    }
                    if ( !ForesterUtil.isEmpty( suffix ) ) {
                        my_name_for_file += suffix;
                    }
                }
                else {
                    if ( !ForesterUtil.isEmpty( default_name ) ) {
                        my_name = new String( default_name );
                    }
                    my_name_for_file = "";
                    if ( !ForesterUtil.isEmpty( full_path ) ) {
                        my_name_for_file = new String( full_path );
                    }
                    else if ( !ForesterUtil.isEmpty( default_name ) ) {
                        my_name_for_file = new String( default_name );
                    }
                    if ( ForesterUtil.isEmpty( my_name_for_file ) ) {
                        if ( !ForesterUtil.isEmpty( phy.getName() ) ) {
                            my_name_for_file = new String( phy.getName() ).replaceAll( " ", "_" );
                        }
                        else if ( phy.getIdentifier() != null ) {
                            final StringBuffer sb = new StringBuffer();
                            if ( !ForesterUtil.isEmpty( phy.getIdentifier().getProvider() ) ) {
                                sb.append( phy.getIdentifier().getProvider() );
                                sb.append( "_" );
                            }
                            sb.append( phy.getIdentifier().getValue() );
                            my_name_for_file = new String( sb.toString().replaceAll( " ", "_" ) );
                        }
                    }
                }
                main_panel.addPhylogenyInNewTab( phy, configuration, my_name, full_path );
                main_panel.getCurrentTreePanel().setTreeFile( new File( my_name_for_file ) );
                lookAtSomeTreePropertiesForAptxControlSettings( phy, main_panel.getControlPanel(), configuration );
            }
            ++i;
        }
    }

    final static void addPhylogenyToPanel( final Phylogeny[] phys,
                                           final Configuration configuration,
                                           final MainPanel main_panel ) {
        final Phylogeny phy = phys[ 0 ];
        main_panel.addPhylogenyInPanel( phy, configuration );
        lookAtSomeTreePropertiesForAptxControlSettings( phy, main_panel.getControlPanel(), configuration );
    }

    final static Color calculateColorFromString( final String str ) {
        final String species_uc = str.toUpperCase();
        char first = species_uc.charAt( 0 );
        char second = ' ';
        char third = ' ';
        if ( species_uc.length() > 1 ) {
            second = species_uc.charAt( 1 );
            if ( species_uc.length() > 2 ) {
                if ( species_uc.indexOf( " " ) > 0 ) {
                    third = species_uc.charAt( species_uc.indexOf( " " ) + 1 );
                }
                else {
                    third = species_uc.charAt( 2 );
                }
            }
        }
        first = Util.normalizeCharForRGB( first );
        second = Util.normalizeCharForRGB( second );
        third = Util.normalizeCharForRGB( third );
        if ( ( first > 235 ) && ( second > 235 ) && ( third > 235 ) ) {
            first = 0;
        }
        else if ( ( first < 80 ) && ( second < 80 ) && ( third < 80 ) ) {
            second = 255;
        }
        return new Color( first, second, third );
    }

    // Returns true if the specified format name can be written
    final static boolean canWriteFormat( final String format_name ) {
        final Iterator<ImageWriter> iter = ImageIO.getImageWritersByFormatName( format_name );
        return iter.hasNext();
    }

    final public static void collapseSpeciesSpecificSubtrees( final Phylogeny phy ) {
        boolean inferred = false;
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.isCollapse() && ( n.getNumberOfDescendants() > 1 ) ) {
                if ( PhylogenyMethods.calculateSumOfDistinctTaxonomies( n ) == 1 ) {
                    Util.collapseSubtree( n, true );
                    if ( !n.getNodeData().isHasTaxonomy() ) {
                        n.getNodeData().setTaxonomy( ( Taxonomy ) n.getAllExternalDescendants().get( 0 ).getNodeData()
                                .getTaxonomy().copy() );
                    }
                    inferred = true;
                }
                else {
                    n.setCollapse( false );
                }
            }
        }
        if ( inferred ) {
            phy.setRerootable( false );
        }
    }

    final static void collapseSubtree( final PhylogenyNode node, final boolean collapse ) {
        node.setCollapse( collapse );
        if ( node.isExternal() ) {
            return;
        }
        final PhylogenyNodeIterator it = new PreorderTreeIterator( node );
        while ( it.hasNext() ) {
            it.next().setCollapse( collapse );
        }
    }

    final static void colorPhylogenyAccordingToConfidenceValues( final Phylogeny tree, final TreePanel tree_panel ) {
        double max_conf = 0.0;
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.getBranchData().setBranchColor( null );
            if ( n.getBranchData().isHasConfidences() ) {
                final double conf = PhylogenyMethods.getConfidenceValue( n );
                if ( conf > max_conf ) {
                    max_conf = conf;
                }
            }
        }
        if ( max_conf > 0.0 ) {
            final Color bg = tree_panel.getTreeColorSet().getBackgroundColor();
            final Color br = tree_panel.getTreeColorSet().getBranchColor();
            for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.getBranchData().isHasConfidences() ) {
                    final double conf = PhylogenyMethods.getConfidenceValue( n );
                    final BranchColor c = new BranchColor( ForesterUtil.calcColor( conf, 0.0, max_conf, bg, br ) );
                    n.getBranchData().setBranchColor( c );
                    final Set<PhylogenyNode> descs = PhylogenyMethods.getAllDescendants( n );
                    for( final PhylogenyNode desc : descs ) {
                        desc.getBranchData().setBranchColor( c );
                    }
                }
            }
        }
    }

    final static void colorPhylogenyAccordingToExternalTaxonomy( final Phylogeny tree, final TreePanel tree_panel ) {
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            it.next().getBranchData().setBranchColor( null );
        }
        for( final PhylogenyNodeIterator it = tree.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.getBranchData().isHasBranchColor() ) {
                final Taxonomy tax = PhylogenyMethods.getExternalDescendantsTaxonomy( n );
                if ( tax != null ) {
                    n.getBranchData().setBranchColor( new BranchColor( tree_panel.calculateTaxonomyBasedColor( tax ) ) );
                    final Set<PhylogenyNode> descs = PhylogenyMethods.getAllDescendants( n );
                    for( final PhylogenyNode desc : descs ) {
                        desc.getBranchData().setBranchColor( new BranchColor( tree_panel
                                .calculateTaxonomyBasedColor( tax ) ) );
                    }
                }
            }
        }
    }

    final static String createDescriptionForTab( final Phylogeny phy, final String full_path ) {
        final StringBuilder desc = new StringBuilder();
        if ( ( phy != null ) && !phy.isEmpty() ) {
            if ( !ForesterUtil.isEmpty( phy.getName() ) ) {
                desc.append( "name: " );
                desc.append( phy.getName() );
                desc.append( "; " );
            }
            if ( phy.getIdentifier() != null ) {
                desc.append( "id: " );
                desc.append( phy.getIdentifier() );
                desc.append( "; " );
            }
            desc.append( "rooted: " );
            desc.append( phy.isRooted() );
            desc.append( "; " );
            desc.append( "rerootable: " );
            desc.append( phy.isRerootable() );
            desc.append( "; " );
            desc.append( "external nodes: " );
            desc.append( phy.getNumberOfExternalNodes() );
            desc.append( "; " );
            desc.append( "branches: " );
            desc.append( phy.getNumberOfBranches() );
            desc.append( "; " );
            desc.append( "maximum descendants: " );
            desc.append( PhylogenyMethods.calculateMaximumNumberOfDescendantsPerNode( phy ) );
            if ( !ForesterUtil.isEmpty( full_path ) && ( full_path.length() <= 50 ) ) {
                desc.append( "; " );
                desc.append( "path: " );
                desc.append( full_path );
            }
        }
        return desc.toString();
    }

    /**
     * Exits with -1.
     * 
     * 
     * @param message
     *            to message to be printed
     */
    final static void dieWithSystemError( final String message ) {
        System.out.println();
        System.out.println( Constants.PRG_NAME + " encountered the following system error: " + message );
        System.out.println( "Please contact the authors." );
        System.out.println( Constants.PRG_NAME + " needs to close." );
        System.out.println();
        System.exit( -1 );
    }

    final static String[] getAvailableFontFamiliesSorted() {
        return AVAILABLE_FONT_FAMILIES_SORTED;
    }

    final static void inferCommonPartOfScientificNames( final Phylogeny tree ) {
        boolean inferred = false;
        for( final PhylogenyNodeIterator it = tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.getNodeData().isHasTaxonomy() && !n.isExternal() ) {
                final String sn = PhylogenyMethods.inferCommonPartOfScientificNameOfDescendants( n );
                if ( !ForesterUtil.isEmpty( sn ) ) {
                    n.getNodeData().setTaxonomy( new Taxonomy() );
                    n.getNodeData().getTaxonomy().setScientificName( sn );
                    inferred = true;
                }
            }
        }
        if ( inferred ) {
            tree.setRerootable( false );
        }
    }

    final static boolean isHasAssignedEvent( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasEvent() ) {
            return false;
        }
        if ( ( node.getNodeData().getEvent() ).isUnassigned() ) {
            return false;
        }
        return true;
    }

    final static boolean isMac() {
        try {
            final String s = ForesterUtil.OS_NAME.toLowerCase();
            return s.startsWith( "mac" );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "minor error: " + e );
            return false;
        }
    }

    final static boolean isJava15() {
        try {
            final String s = ForesterUtil.JAVA_VERSION;
            return s.startsWith( "1.5" );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "minor error: " + e );
            return false;
        }
    }

    final static boolean isUsOrCanada() {
        try {
            if ( ( Locale.getDefault().equals( Locale.CANADA ) ) || ( Locale.getDefault().equals( Locale.US ) ) ) {
                return true;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return false;
    }

    final static boolean isWindows() {
        try {
            final String s = ForesterUtil.OS_NAME.toLowerCase();
            return s.indexOf( "win" ) > -1;
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "minor error: " + e );
            return false;
        }
    }

    final static void launchWebBrowser( final URI uri,
                                        final boolean is_applet,
                                        final JApplet applet,
                                        final String frame_name ) throws IOException {
        if ( is_applet ) {
            applet.getAppletContext().showDocument( uri.toURL(), frame_name );
        }
        else {
            // This requires Java 1.6:
            // =======================
            // boolean no_desktop = false;
            // try {
            // if ( Desktop.isDesktopSupported() ) {
            // System.out.println( "desktop supported" );
            // final Desktop dt = Desktop.getDesktop();
            // dt.browse( uri );
            // }
            // else {
            // no_desktop = true;
            // }
            // }
            // catch ( final Exception ex ) {
            // ex.printStackTrace();
            // no_desktop = true;
            // }
            // catch ( final Error er ) {
            // er.printStackTrace();
            // no_desktop = true;
            // }
            // if ( no_desktop ) {
            // System.out.println( "desktop not supported" );
            try {
                openUrlInWebBrowser( uri.toString() );
            }
            catch ( final Exception e ) {
                throw new IOException( e );
            }
            // }
        }
    }

    final static void lookAtSomeTreePropertiesForAptxControlSettings( final Phylogeny t,
                                                                      final ControlPanel atv_control,
                                                                      final Configuration configuration ) {
        if ( ( t != null ) && !t.isEmpty() ) {
            if ( !ForesterUtil.isHasAtLeastOneBranchLengthLargerThanZero( t ) ) {
                atv_control.setDrawPhylogram( false );
                atv_control.setDrawPhylogramEnabled( false );
            }
            if ( configuration.doGuessCheckOption( Configuration.display_as_phylogram ) ) {
                if ( atv_control.getDisplayAsPhylogramCb() != null ) {
                    if ( ForesterUtil.isHasAtLeastOneBranchLengthLargerThanZero( t ) ) {
                        atv_control.setDrawPhylogram( true );
                        atv_control.setDrawPhylogramEnabled( true );
                    }
                    else {
                        atv_control.setDrawPhylogram( false );
                    }
                }
            }
            if ( configuration.doGuessCheckOption( Configuration.write_confidence_values ) ) {
                if ( atv_control.getWriteConfidenceCb() != null ) {
                    if ( ForesterUtil.isHasAtLeastOneBranchWithSupportValues( t ) ) {
                        atv_control.setCheckbox( Configuration.write_confidence_values, true );
                    }
                    else {
                        atv_control.setCheckbox( Configuration.write_confidence_values, false );
                    }
                }
            }
            if ( configuration.doGuessCheckOption( Configuration.write_events ) ) {
                if ( atv_control.getShowEventsCb() != null ) {
                    if ( ForesterUtil.isHasAtLeastNodeWithEvent( t ) ) {
                        atv_control.setCheckbox( Configuration.write_events, true );
                    }
                    else {
                        atv_control.setCheckbox( Configuration.write_events, false );
                    }
                }
            }
        }
    }

    final private static char normalizeCharForRGB( char c ) {
        c -= 65;
        c *= 10.2;
        c = c > 255 ? 255 : c;
        c = c < 0 ? 0 : c;
        return c;
    }

    final private static void openUrlInWebBrowser( final String url ) throws IOException, ClassNotFoundException,
            SecurityException, NoSuchMethodException, IllegalArgumentException, IllegalAccessException,
            InvocationTargetException, InterruptedException {
        final String os = System.getProperty( "os.name" );
        final Runtime runtime = Runtime.getRuntime();
        if ( os.toLowerCase().startsWith( "win" ) ) {
            Runtime.getRuntime().exec( "rundll32 url.dll,FileProtocolHandler " + url );
        }
        else if ( isMac() ) {
            final Class<?> file_mgr = Class.forName( "com.apple.eio.FileManager" );
            final Method open_url = file_mgr.getDeclaredMethod( "openURL", new Class[] { String.class } );
            open_url.invoke( null, new Object[] { url } );
        }
        else {
            final String[] browsers = { "firefox", "opera", "konqueror", "mozilla", "netscape", "epiphany" };
            String browser = null;
            for( int i = 0; ( i < browsers.length ) && ( browser == null ); ++i ) {
                if ( runtime.exec( new String[] { "which", browsers[ i ] } ).waitFor() == 0 ) {
                    browser = browsers[ i ];
                }
            }
            if ( browser == null ) {
                throw new IOException( "could not find a web browser to open [" + url + "] in" );
            }
            else {
                runtime.exec( new String[] { browser, url } );
            }
        }
    }

    final static void openWebsite( final String url, final boolean is_applet, final JApplet applet ) throws IOException {
        try {
            Util.launchWebBrowser( new URI( url ), is_applet, applet, Constants.PRG_NAME );
        }
        catch ( final Exception e ) {
            throw new IOException( e );
        }
    }

    final static void printAppletMessage( final String applet_name, final String message ) {
        System.out.println( "[" + applet_name + "] > " + message );
    }

    final static Phylogeny[] readPhylogenies( final PhylogenyParser parser, final File file ) throws IOException {
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final Phylogeny[] trees = factory.create( file, parser );
        if ( ( trees == null ) || ( trees.length == 0 ) ) {
            throw new PhylogenyParserException( "Unable to parse phylogeny from file: " + file );
        }
        return trees;
    }

    final static Phylogeny[] readPhylogeniesFromUrl( final URL url, final boolean phyloxml_validate_against_xsd )
            throws FileNotFoundException, IOException {
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        PhylogenyParser parser = null;
        if ( url.getHost().toLowerCase().indexOf( "tolweb" ) >= 0 ) {
            parser = new TolParser();
        }
        else {
            parser = ForesterUtil.createParserDependingOnUrlContents( url, phyloxml_validate_against_xsd );
        }
        return factory.create( url.openStream(), parser );
    }

    final static void removeBranchColors( final Phylogeny phy ) {
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            it.next().getBranchData().setBranchColor( null );
        }
    }

    final static void showErrorMessage( final Component parent, final String error_msg ) {
        printAppletMessage( Constants.PRG_NAME, error_msg );
        JOptionPane.showMessageDialog( parent, error_msg, "[" + Constants.PRG_NAME + " " + Constants.VERSION
                + "] Error", JOptionPane.ERROR_MESSAGE );
    }

    final static void unexpectedError( final Error err ) {
        err.printStackTrace();
        final StringBuffer sb = new StringBuffer();
        for( final StackTraceElement s : err.getStackTrace() ) {
            sb.append( s + "\n" );
        }
        JOptionPane
                .showMessageDialog( null,
                                    "An unexpected (possibly severe) error has occured - terminating. \nPlease contact: "
                                            + Constants.AUTHOR_EMAIL + " \nError: " + err + "\n" + sb,
                                    "Unexpected Severe Error [" + Constants.PRG_NAME + " " + Constants.VERSION + "]",
                                    JOptionPane.ERROR_MESSAGE );
        System.exit( -1 );
    }

    final static void unexpectedException( final Exception ex ) {
        ex.printStackTrace();
        final StringBuffer sb = new StringBuffer();
        for( final StackTraceElement s : ex.getStackTrace() ) {
            sb.append( s + "\n" );
        }
        JOptionPane.showMessageDialog( null, "An unexpected exception has occured. \nPlease contact: "
                + Constants.AUTHOR_EMAIL + " \nException: " + ex + "\n" + sb, "Unexpected Exception ["
                + Constants.PRG_NAME + Constants.VERSION + "]", JOptionPane.ERROR_MESSAGE );
    }

    final static String writePhylogenyToGraphicsFile( final String file_name,
                                                      int width,
                                                      int height,
                                                      final TreePanel tree_panel,
                                                      final ControlPanel ac,
                                                      final GraphicsExportType type,
                                                      final Options options ) throws IOException {
        if ( !options.isGraphicsExportUsingActualSize() ) {
            if ( options.isGraphicsExportVisibleOnly() ) {
                throw new IllegalArgumentException( "cannot export visible rectangle only without exporting in actual size" );
            }
            tree_panel.setParametersForPainting( options.getPrintSizeX(), options.getPrintSizeY(), true );
            tree_panel.resetPreferredSize();
            tree_panel.repaint();
        }
        final RenderingHints rendering_hints = new RenderingHints( RenderingHints.KEY_RENDERING,
                                                                   RenderingHints.VALUE_RENDER_QUALITY );
        rendering_hints.put( RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY );
        if ( options.isAntialiasPrint() ) {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        }
        else {
            rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
            rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF );
        }
        final Phylogeny phylogeny = tree_panel.getPhylogeny();
        if ( ( phylogeny == null ) || phylogeny.isEmpty() ) {
            return "";
        }
        final File file = new File( file_name );
        if ( file.isDirectory() ) {
            throw new IOException( "\"" + file_name + "\" is a directory" );
        }
        Rectangle visible = null;
        if ( !options.isGraphicsExportUsingActualSize() ) {
            width = options.getPrintSizeX();
            height = options.getPrintSizeY();
        }
        else if ( options.isGraphicsExportVisibleOnly() ) {
            visible = tree_panel.getVisibleRect();
            width = visible.width;
            height = visible.height;
        }
        final BufferedImage buffered_img = new BufferedImage( width, height, BufferedImage.TYPE_INT_RGB );
        Graphics2D g2d = buffered_img.createGraphics();
        g2d.setRenderingHints( rendering_hints );
        int x = 0;
        int y = 0;
        if ( options.isGraphicsExportVisibleOnly() ) {
            g2d = ( Graphics2D ) g2d.create( -visible.x, -visible.y, visible.width, visible.height );
            g2d.setClip( null );
            x = visible.x;
            y = visible.y;
        }
        tree_panel.paintPhylogeny( g2d, false, true, width, height, x, y );
        if ( type == GraphicsExportType.TIFF ) {
            writeToTiff( file, buffered_img );
        }
        else {
            ImageIO.write( buffered_img, type.toString(), file );
        }
        g2d.dispose();
        System.gc();
        if ( !options.isGraphicsExportUsingActualSize() ) {
            tree_panel.getMainPanel().getControlPanel().showWhole();
        }
        String msg = file.toString();
        if ( ( width > 0 ) && ( height > 0 ) ) {
            msg += " [size: " + width + ", " + height + "]";
        }
        return msg;
    }

    final static void writeToTiff( final File file, final BufferedImage image ) throws IOException {
        // See: http://log.robmeek.com/2005/08/write-tiff-in-java.html
        ImageWriter writer = null;
        ImageOutputStream ios = null;
        // Find an appropriate writer:
        final Iterator<ImageWriter> it = ImageIO.getImageWritersByFormatName( "TIF" );
        if ( it.hasNext() ) {
            writer = it.next();
        }
        else {
            throw new IOException( "failed to get TIFF image writer" );
        }
        // Setup writer:
        ios = ImageIO.createImageOutputStream( file );
        writer.setOutput( ios );
        final ImageWriteParam image_write_param = new ImageWriteParam( Locale.getDefault() );
        image_write_param.setCompressionMode( ImageWriteParam.MODE_EXPLICIT );
        // see writeParam.getCompressionTypes() for available compression type
        // strings.
        image_write_param.setCompressionType( "PackBits" );
        final String t[] = image_write_param.getCompressionTypes();
        for( final String string : t ) {
            System.out.println( string );
        }
        // Convert to an IIOImage:
        final IIOImage iio_image = new IIOImage( image, null, null );
        writer.write( null, iio_image, image_write_param );
    }

    // See: http://www.xml.nig.ac.jp/tutorial/rest/index.html#2.2
    // static void openDDBJRest() throws IOException {
    // //set URL
    // URL url = new URL( "http://xml.nig.ac.jp/rest/Invoke" );
    // //set parameter
    // String query = "service=GetEntry&method=getDDBJEntry&accession=AB000100";
    // //make connection
    // URLConnection urlc = url.openConnection();
    // //use post mode
    // urlc.setDoOutput( true );
    // urlc.setAllowUserInteraction( false );
    // //send query
    // PrintStream ps = new PrintStream( urlc.getOutputStream() );
    // ps.print( query );
    // ps.close();
    // //get result
    // BufferedReader br = new BufferedReader( new InputStreamReader(
    // urlc.getInputStream() ) );
    // String l = null;
    // while ( ( l = br.readLine() ) != null ) {
    // System.out.println( l );
    // }
    // br.close();
    // }
    static enum GraphicsExportType {
        GIF( "gif" ), JPG( "jpg" ), PDF( "pdf" ), PNG( "png" ), TIFF( "tif" ), BMP( "bmp" );

        private final String _suffix;

        private GraphicsExportType( final String suffix ) {
            _suffix = suffix;
        }

        @Override
        public String toString() {
            return _suffix;
        }
    }
}
