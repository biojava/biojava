// $Id: rio.java,v 1.11 2010/01/16 02:15:34 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.application;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.sdi.DistanceCalculator;
import org.forester.sdi.RIO;
import org.forester.sdi.SDIR;
import org.forester.util.ForesterUtil;

public class rio {

    final static private String  PRG_NAME                         = "RIO";
    final static private String  PRG_VERSION                      = "2.03 ALPHA";
    final static private String  PRG_DATE                         = "2010.01.15";
    final static private String  E_MAIL                           = "czmasek@burnham.org";
    final static private String  WWW                              = "www.phylosoft.org/forester/";
    final static private boolean TIME                             = true;
    final static private boolean VERBOSE                          = true;
    // For method getDistances -- calculation of distances.
    final static private boolean MINIMIZE_COST                    = false;
    // For method getDistances -- calculation of distances.
    final static private boolean MINIMIZE_DUPS                    = true;
    // For method getDistances -- calculation of distances.
    final static private boolean MINIMIZE_HEIGHT                  = true;
    final static private int     WARN_NO_ORTHOS_DEFAULT           = 2;
    final static private int
                                 // How many sd away from mean to root.
                                 WARN_MORE_THAN_ONE_ORTHO_DEFAULT = 2;
    // How many sd away from mean to LCA of orthos.
    final static private double  THRESHOLD_ULTRA_PARALOGS_DEFAULT = 50;
    // How many sd away from mean to LCA of orthos.
    final static private double  WARN_ONE_ORTHO_DEFAULT           = 2;

    // Factor between the two distances to their LCA
    // (larger/smaller).
    // Factor between the two distances to their LCA
    // (larger/smaller).
    /**
     * Calculates the mean and standard deviation of all nodes of Phylogeny t
     * which have a bootstrap values zero or more. Returns null in case of
     * failure (e.g t has no bootstrap values, or just one).
     * <p>
     * 
     * @param t
     *            reference to a tree with bootstrap values
     * @return Array of doubles, [0] is the mean, [1] the standard deviation
     */
    private static double[] calculateMeanBoostrapValue( final Phylogeny t ) {
        double b = 0;
        int n = 0;
        long sum = 0;
        double x = 0.0, mean = 0.0;
        final double[] da = new double[ 2 ];
        final Vector<Double> bv = new Vector<Double>();
        PhylogenyNode node = null;
        PreorderTreeIterator i = null;
        i = new PreorderTreeIterator( t );
        // Calculates the mean.
        while ( i.hasNext() ) {
            node = i.next();
            if ( !( ( node.getParent() != null ) && node.getParent().isRoot()
                    && ( PhylogenyMethods.getConfidenceValue( node.getParent().getChildNode1() ) > 0 )
                    && ( PhylogenyMethods.getConfidenceValue( node.getParent().getChildNode2() ) > 0 ) && ( node
                    .getParent().getChildNode2() == node ) ) ) {
                b = PhylogenyMethods.getConfidenceValue( node );
                if ( b > 0 ) {
                    sum += b;
                    bv.addElement( new Double( b ) );
                    n++;
                }
            }
            // i.next();
        }
        if ( n < 2 ) {
            return null;
        }
        mean = ( double ) sum / n;
        // Calculates the standard deviation.
        sum = 0;
        for( int j = 0; j < n; ++j ) {
            b = ( bv.elementAt( j ) ).intValue();
            x = b - mean;
            sum += ( x * x );
        }
        da[ 0 ] = mean;
        da[ 1 ] = java.lang.Math.sqrt( sum / ( n - 1.0 ) );
        return da;
    }

    private final static void errorInCommandLine() {
        System.out.println( "\nrio: Error in command line.\n" );
        printHelp();
        System.exit( -1 );
    }

    // Uses DistanceCalculator to calculate distances.
    private final static StringBuffer getDistances( final File tree_file_for_dist_val,
                                                    final File outfile,
                                                    final Phylogeny species_tree,
                                                    final String seq_name,
                                                    final ArrayList<String> al_ortholog_names_for_dc,
                                                    final HashMap<String, Integer> ortholog_hashmap,
                                                    final HashMap<String, Integer> super_ortholog_hashmap,
                                                    final int warn_more_than_one_ortho,
                                                    final int warn_no_orthos,
                                                    final double warn_one_ortho,
                                                    final int bootstraps,
                                                    final double t_orthologs_dc ) throws IOException {
        Phylogeny consensus_tree = null;
        Phylogeny
        // to be a consensus tree.
        assigned_cons_tree = null;
        final SDIR sdiunrooted = new SDIR();
        final ArrayList<PhylogenyNode> al_ortholog_nodes = new ArrayList<PhylogenyNode>();
        double m = 0.0;
        double sd = 0.0;
        double d = 0.0;
        int n = 0;
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        consensus_tree = factory.create( tree_file_for_dist_val, new PhyloXmlParser() )[ 0 ];
        PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, consensus_tree );
        assigned_cons_tree = sdiunrooted.infer( consensus_tree,
                                                species_tree,
                                                rio.MINIMIZE_COST,
                                                rio.MINIMIZE_DUPS,
                                                rio.MINIMIZE_HEIGHT,
                                                true,
                                                1 )[ 0 ];
        final DistanceCalculator dc = new DistanceCalculator();
        final StringBuffer sb = new StringBuffer();
        sb.append( "Given the threshold for distance calculations (" + ForesterUtil.roundToInt( t_orthologs_dc )
                + "): " );
        // No orthologs.
        if ( al_ortholog_names_for_dc.size() == 0 ) {
            dc.setTree( assigned_cons_tree );
            // Remark. Calculation of mean and sd _does_ include the node
            // with seq_name.
            m = dc.getMean();
            sd = dc.getStandardDeviation();
            d = dc.getDistanceToRoot( seq_name );
            n = dc.getN();
            sb.append( "No sequence is considered orthologous to query."
                    + "\ndistance of query to root                     = " + ForesterUtil.FORMATTER_06.format( d )
                    + "\nmean of distances (for all sequences) to root = " + ForesterUtil.FORMATTER_06.format( m )
                    + "\nsd of distances (for all sequences) to root   = " + ForesterUtil.FORMATTER_06.format( sd )
                    + "\nn (sum of sequences in alignment plus query)  = " + n );
            if ( !( ( ( m - ( warn_no_orthos * sd ) ) < d ) && ( ( m + ( warn_no_orthos * sd ) ) > d ) ) ) {
                sb.append( "\nWARNING: distance of query to root is outside of mean+/-" + warn_no_orthos + "*sd!" );
            }
        }
        // One ortholog.
        else if ( al_ortholog_names_for_dc.size() == 1 ) {
            final String name_of_ortholog = al_ortholog_names_for_dc.get( 0 );
            al_ortholog_nodes.add( assigned_cons_tree.getNode( name_of_ortholog ) );
            al_ortholog_nodes.add( assigned_cons_tree.getNode( seq_name ) );
            dc.setTreeAndExtNodes( assigned_cons_tree, al_ortholog_nodes );
            // Remark. Calculation of mean _does_ include the node
            // with seq_name.
            d = dc.getDistanceToLCA( seq_name );
            final double d_o = dc.getDistanceToLCA( name_of_ortholog );
            sb.append( "One sequence is considered orthologous to query." + "\nLCA is LCA of query and its ortholog."
                    + "\ndistance of query to LCA    = " + ForesterUtil.FORMATTER_06.format( d )
                    + "\ndistance of ortholog to LCA = " + ForesterUtil.FORMATTER_06.format( d_o ) );
            if ( ( d_o > 0.0 )
                    && ( d > 0.0 )
                    && ( ( ( d_o >= d ) && ( ( d_o / d ) > warn_one_ortho ) ) || ( ( d_o < d ) && ( ( d / d_o ) > warn_one_ortho ) ) ) ) {
                sb.append( "\nWARNING: Ratio of distances to LCA is greater than " + warn_one_ortho + "!" );
            }
            else if ( ( ( d_o == 0.0 ) || ( d == 0.0 ) ) && ( ( d_o != 0.0 ) || ( d != 0.0 ) ) ) {
                sb.append( "\nWARNING: Ratio could not be calculated, " + " one distance is 0.0!" );
            }
        }
        // More than one ortholog.
        else {
            for( int i = 0; i < al_ortholog_names_for_dc.size(); ++i ) {
                al_ortholog_nodes.add( assigned_cons_tree.getNodeViaSequenceName( al_ortholog_names_for_dc.get( i ) ) );
            }
            al_ortholog_nodes.add( assigned_cons_tree.getNodesViaSequenceName( seq_name ).get( 0 ) );
            dc.setTreeAndExtNodes( assigned_cons_tree, al_ortholog_nodes );
            // Remark. Calculation of mean and sd _does_ include the node
            // with seq_name.
            m = dc.getMean();
            sd = dc.getStandardDeviation();
            d = dc.getDistanceToLCA( seq_name );
            n = dc.getN();
            sb.append( "More than one sequence is considered orthologous to query."
                    + "\nLCA is LCA of query and its orthologs."
                    + "\ndistance of query to LCA                               = "
                    + ForesterUtil.FORMATTER_06.format( d )
                    + "\nmean of distances (for query and its orthologs) to LCA = "
                    + ForesterUtil.FORMATTER_06.format( m )
                    + "\nsd of distances (for query and its orthologs) to LCA   = "
                    + ForesterUtil.FORMATTER_06.format( sd )
                    + "\nn (sum of orthologs plus query)                        = " + n );
            if ( !( ( ( m - ( warn_more_than_one_ortho * sd ) ) < d ) && ( ( m + ( warn_more_than_one_ortho * sd ) ) > d ) ) ) {
                sb.append( "\n!WARNING: distance of query to LCA is outside of mean+/-" + warn_more_than_one_ortho
                        + "*sd!" );
            }
        }
        return sb;
    }

    public static void main( final String[] args ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        File species_tree_file = null;
        File multiple_trees_file = null;
        File outfile = null;
        File distance_matrix_file = null;
        File tree_file_for_dist_val = null;
        File tree_file_for_avg_bs = null;
        String seq_name = "";
        String arg = "";
        boolean output_ultraparalogs = false;
        ArrayList<String> orthologs_al_for_dc = null;
        double t_orthologs = 0.0;
        double t_sn = 0.0;
        double t_orthologs_dc = 0.0;
        double[] bs_mean_sd = null;
        int sort = 13;
        Phylogeny species_tree = null;
        RIO rio_instance = null;
        PrintWriter out = null;
        long time = 0;
        int warn_no_orthos = WARN_NO_ORTHOS_DEFAULT;
        int warn_more_than_one_ortho = WARN_MORE_THAN_ONE_ORTHO_DEFAULT;
        double warn_one_ortho = WARN_ONE_ORTHO_DEFAULT;
        double threshold_ultra_paralogs = THRESHOLD_ULTRA_PARALOGS_DEFAULT;
        if ( args.length < 2 ) {
            printHelp();
            System.exit( 0 );
        }
        else if ( ( args.length < 3 ) || ( args.length > 18 ) ) {
            errorInCommandLine();
        }
        for( int i = 0; i < args.length; ++i ) {
            if ( args[ i ].trim().charAt( 0 ) != 'p' ) {
                if ( args[ i ].trim().length() < 3 ) {
                    errorInCommandLine();
                }
                else {
                    arg = args[ i ].trim().substring( 2 );
                }
            }
            try {
                switch ( args[ i ].trim().charAt( 0 ) ) {
                    case 'M':
                        multiple_trees_file = new File( arg );
                        break;
                    case 'N':
                        seq_name = arg;
                        break;
                    case 'S':
                        species_tree_file = new File( arg );
                        break;
                    case 'O':
                        outfile = new File( arg );
                        break;
                    case 'D':
                        distance_matrix_file = new File( arg );
                        break;
                    case 'T':
                        tree_file_for_dist_val = new File( arg );
                        break;
                    case 't':
                        tree_file_for_avg_bs = new File( arg );
                        break;
                    case 'p':
                        output_ultraparalogs = true;
                        break;
                    case 'P':
                        sort = Integer.parseInt( arg );
                        if ( ( sort < 0 ) || ( sort > 17 ) ) {
                            errorInCommandLine();
                        }
                        break;
                    case 'L':
                        t_orthologs = Double.parseDouble( arg );
                        break;
                    case 'B':
                        t_sn = Double.parseDouble( arg );
                        break;
                    case 'U':
                        t_orthologs_dc = Double.parseDouble( arg );
                        break;
                    case 'v':
                        threshold_ultra_paralogs = Double.parseDouble( arg );
                        break;
                    case 'X':
                        warn_more_than_one_ortho = Integer.parseInt( arg );
                        break;
                    case 'Y':
                        warn_no_orthos = Integer.parseInt( arg );
                        break;
                    case 'Z':
                        warn_one_ortho = Double.parseDouble( arg );
                        break;
                    default:
                        errorInCommandLine();
                }
            }
            catch ( final Exception e ) {
                errorInCommandLine();
            }
        }
        if ( ( seq_name == "" ) || ( species_tree_file == null ) || ( multiple_trees_file == null )
                || ( outfile == null ) ) {
            errorInCommandLine();
        }
        if ( ( sort < 0 ) || ( sort > 17 ) ) {
            errorInCommandLine();
        }
        if ( ( sort > 2 ) && ( distance_matrix_file == null ) ) {
            errorInCommandLine();
        }
        if ( VERBOSE ) {
            System.out.println( "\nMultiple trees file:                          " + multiple_trees_file );
            System.out.println( "Seq name:                                     " + seq_name );
            System.out.println( "Species tree file:                            " + species_tree_file );
            System.out.println( "Outfile:                                      " + outfile );
            if ( distance_matrix_file != null ) {
                System.out.println( "Distance matrix file:                         " + distance_matrix_file );
            }
            if ( tree_file_for_dist_val != null ) {
                if ( tree_file_for_avg_bs == null ) {
                    System.out.println( "Phy to read dists and calc mean support from: " + tree_file_for_dist_val );
                }
                else {
                    System.out.println( "Phylogeny to read dist values from:                " + tree_file_for_dist_val );
                }
            }
            if ( tree_file_for_avg_bs != null ) {
                System.out.println( "Phylogeny to calc mean bootstrap from:             " + tree_file_for_avg_bs );
            }
            System.out.println( "Sort:                                         " + sort );
            System.out.println( "Threshold orthologs:                          " + t_orthologs );
            System.out.println( "Threshold subtree neighborings:               " + t_sn );
            System.out.println( "Threshold orthologs for distance calc.:       " + t_orthologs_dc );
            if ( output_ultraparalogs ) {
                System.out.println( "Threshold ultra paralogs:                     " + threshold_ultra_paralogs );
            }
            System.out.println( "More than one ortholog  sd diff:              " + warn_more_than_one_ortho );
            System.out.println( "No  orthologs           sd diff:              " + warn_no_orthos );
            System.out.println( "One ortholog            factor :              " + warn_one_ortho + "\n" );
        }
        if ( TIME && VERBOSE ) {
            time = System.currentTimeMillis();
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            species_tree = factory.create( species_tree_file, new PhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        if ( !species_tree.isRooted() ) {
            ForesterUtil.printErrorMessage( PRG_NAME, "Species tree is not rooted" );
            System.exit( -1 );
        }
        if ( !species_tree.isCompletelyBinary() ) {
            ForesterUtil.printErrorMessage( PRG_NAME, "Species tree is not completely binary" );
            System.exit( -1 );
        }
        rio_instance = new RIO();
        final StringBuffer output = new StringBuffer();
        try {
            if ( distance_matrix_file != null ) {
                rio_instance.readDistanceMatrix( distance_matrix_file );
            }
            rio_instance.inferOrthologs( multiple_trees_file, species_tree.copy(), seq_name );
            output.append( rio_instance.inferredOrthologsToString( seq_name, sort, t_orthologs, t_sn ) );
            if ( tree_file_for_dist_val != null ) {
                orthologs_al_for_dc = rio_instance.inferredOrthologsToArrayList( seq_name, t_orthologs_dc );
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                if ( tree_file_for_avg_bs != null ) {
                    final Phylogeny p = factory.create( tree_file_for_avg_bs, new PhyloXmlParser() )[ 0 ];
                    bs_mean_sd = calculateMeanBoostrapValue( p );
                }
                else {
                    final Phylogeny p = factory.create( tree_file_for_dist_val, new PhyloXmlParser() )[ 0 ];
                    bs_mean_sd = calculateMeanBoostrapValue( p );
                }
                if ( ( bs_mean_sd != null ) && ( bs_mean_sd.length == 2 ) ) {
                    final double bs_mean = bs_mean_sd[ 0 ];
                    final double bs_sd = bs_mean_sd[ 1 ];
                    output.append( "\n\nMean bootstrap value of consensus tree (sd): "
                            + ForesterUtil.roundToInt( ( bs_mean * 100.0 ) / rio_instance.getBootstraps() ) + "% (+/-"
                            + ForesterUtil.roundToInt( ( bs_sd * 100.0 ) / rio_instance.getBootstraps() ) + "%)\n" );
                }
                output.append( "\n\nDistance values:\n" );
                output.append( getDistances( tree_file_for_dist_val,
                                             outfile,
                                             species_tree,
                                             seq_name,
                                             orthologs_al_for_dc,
                                             rio_instance.getInferredOrthologs( seq_name ),
                                             rio_instance.getInferredSuperOrthologs( seq_name ),
                                             warn_more_than_one_ortho,
                                             warn_no_orthos,
                                             warn_one_ortho,
                                             rio_instance.getBootstraps(),
                                             t_orthologs_dc ) );
            }
            if ( output_ultraparalogs ) {
                output.append( "\n\nUltra paralogs:\n" );
                output.append( rio_instance
                        .inferredUltraParalogsToString( seq_name, sort > 2, threshold_ultra_paralogs ) );
            }
            output.append( "\n\nSort priority: " + RIO.getOrder( sort ) );
            output.append( "\nExt nodes    : " + rio_instance.getExtNodesOfAnalyzedGeneTrees() );
            output.append( "\nSamples      : " + rio_instance.getBootstraps() + "\n" );
            out = new PrintWriter( new FileWriter( outfile ), true );
        }
        catch ( final Exception e ) {
            ForesterUtil.printErrorMessage( PRG_NAME, e.getLocalizedMessage() );
            e.printStackTrace();
            System.exit( -1 );
        }
        out.println( output );
        out.close();
        ForesterUtil.programMessage( PRG_NAME, "wrote results to \"" + outfile + "\"" );
        if ( TIME && VERBOSE ) {
            time = System.currentTimeMillis() - time;
            ForesterUtil.programMessage( PRG_NAME, "time: " + time + "ms" );
        }
        ForesterUtil.programMessage( PRG_NAME, "OK." );
        System.exit( 0 );
    }

    private final static void printHelp() {
        System.out.println( "M= (String) Multiple gene tree file (mandatory)" );
        System.out.println( "N= (String) Query sequence name (mandatory)" );
        System.out.println( "S= (String) Species tree file (mandatory)" );
        System.out.println( "O= (String) Output file name -- overwritten without warning! (mandatory)" );
        System.out.println( "D= (String) Distance matrix file for pairwise distances" );
        System.out.println( "T= (String) Phylogeny file for distances of query to LCA" );
        System.out.println( "            of orthologs and for mean bootstrap value (if t= is not used)," );
        System.out.println( "            must be binary )" );
        System.out.println( "t= (String) Phylogeny file for mean bootstrap value (if this option is used," );
        System.out.println( "            the mean bootstrap value is not calculated from the tree read in" );
        System.out.println( "            with T=), not necessary binary" );
        System.out.println( "p           To output ultra paralogs" );
        System.out.println( "P= (int)    Sort priority" );
        System.out.println( "L= (double) Threshold orthologs for output" );
        System.out.println( "U= (double) Threshold orthologs for distance calculation" );
        System.out.println( "X= (int)    More than one ortholog: " );
        System.out.println( "            numbers of sd the dist. to LCA has to differ from mean to generate a warning" );
        System.out.println( "Y= (int)    No orthologs:" );
        System.out.println( "            Numbers of sd the dist to root has to differ from mean to generate a warning" );
        System.out.println( "Z= (double) One ortholog:" );
        System.out.println( "            threshold for factor between the two distances to their LCA (larger/smaller)" );
        System.out.println( "            to generate a warning" );
        System.out.println();
        System.out.println( " Sort priority (\"P=\"):" );
        System.out.println( RIO.getOrderHelp().toString() );
        System.out.println();
        System.out
                .println( " Example: \"rio M=gene_trees.xml N=bcl2_NEMVE S=species_tree.xml D=distances P=13 p O=out\"" );
        System.out.println();
    }
}
