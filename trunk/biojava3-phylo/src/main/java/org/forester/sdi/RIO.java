// $Id: RIO.java,v 1.26 2010/01/16 02:15:34 cmzmasek Exp $
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

package org.forester.sdi;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogenyinference.DistanceMatrix;
import org.forester.phylogenyinference.SymmetricalDistanceMatrixParser;
import org.forester.util.ForesterUtil;

/*
 * @author Christian M. Zmasek
 */
public final class RIO {

    private final static boolean                      ROOT_BY_MINIMIZING_MAPPING_COST = false;
    private final static boolean                      ROOT_BY_MINIMIZING_SUM_OF_DUPS  = true;
    private final static boolean                      ROOT_BY_MINIMIZING_TREE_HEIGHT  = true;
    private final static boolean                      TIME                            = false;
    private HashMap<String, HashMap<String, Integer>> _o_hash_maps;
    private HashMap<String, HashMap<String, Integer>> _so_hash_maps;
    private HashMap<String, HashMap<String, Integer>> _up_hash_maps;
    private HashMap<String, HashMap<String, Integer>> _sn_hash_maps;                          // HashMap of HashMaps
    private DistanceMatrix                            _m;
    private HashMap<String, Double>                   _l;
    private String[]                                  _seq_names;
    private int                                       _bootstraps;
    private int                                       _ext_nodes_;
    private long                                      _time;

    /**
     * Default constructor.
     */
    public RIO() {
        reset();
    }

    /**
     * Returns the numbers of trees analyzed.
     * 
     * @return the numbers of trees analyzed
     */
    public final int getBootstraps() {
        return _bootstraps;
    }

    // Helper method for inferredOrthologsToString.
    // inferredOrthologsToArrayList,
    // and inferredUltraParalogsToString.
    private final double getBootstrapValueFromHash( final HashMap<String, Integer> h, final String name ) {
        if ( !h.containsKey( name ) ) {
            return 0.0;
        }
        final int i = h.get( name );
        return ( i * 100.0 / getBootstraps() );
    }

    /**
     * Returns the distance to a sequences/taxa after a distance list file has
     * been read in with readDistanceList(File). Throws an exception if name is
     * not found or if no list has been read in.
     * 
     * @param name
     *            a sequence name
     */
    public final double getDistance( String name ) {
        double distance = 0.0;
        name = name.trim();
        if ( _l == null ) {
            throw new IllegalStateException( "Distance list has probably not been read in (successfully)." );
        }
        if ( _l.get( name ) == null ) {
            throw new IllegalArgumentException( name + " not found." );
        }
        distance = ( _l.get( name ) ).doubleValue();
        return distance;
    }

    public final double getDistance( final String name1, final String name2 ) {
        try {
            return _m.getValue( _m.getIndex( name1 ), _m.getIndex( name2 ) );
        }
        catch ( final Exception e ) {
            return 1;
        }
    }

    /**
     * Returns the numbers of number of ext nodes in gene trees analyzed (after
     * stripping).
     * 
     * @return number of ext nodes in gene trees analyzed (after stripping)
     */
    public final int getExtNodesOfAnalyzedGeneTrees() {
        return _ext_nodes_;
    }

    /**
     * Returns a HashMap containing the inferred orthologs of the external gene
     * tree node with the sequence name seq_name. Sequence names are the keys
     * (String), numbers of observations are the values (Int). Orthologs are to
     * be inferred by method "inferOrthologs". Throws an exception if seq_name
     * is not found.
     * 
     * @param seq_name
     *            sequence name of a external node of the gene trees
     * @return HashMap containing the inferred orthologs
     *         (name(String)->value(Int))
     */
    public final HashMap<String, Integer> getInferredOrthologs( final String seq_name ) {
        if ( _o_hash_maps == null ) {
            return null;
        }
        return _o_hash_maps.get( seq_name );
    }

    private final HashMap<String, Integer> getInferredSubtreeNeighbors( final String seq_name ) {
        if ( _sn_hash_maps == null ) {
            return null;
        }
        return _sn_hash_maps.get( seq_name );
    }

    /**
     * Returns a HashMap containing the inferred "super orthologs" of the
     * external gene tree node with the sequence name seq_name. Sequence names
     * are the keys (String), numbers of observations are the values (Int).
     * Super orthologs are to be inferred by method "inferOrthologs". Throws an
     * exception if seq_name is not found.
     * 
     * @param seq_name
     *            sequence name of a external node of the gene trees
     * @return HashMap containing the inferred super orthologs
     *         (name(String)->value(Int))
     */
    public final HashMap<String, Integer> getInferredSuperOrthologs( final String seq_name ) {
        if ( _so_hash_maps == null ) {
            return null;
        }
        return _so_hash_maps.get( seq_name );
    }

    /**
     * Returns a HashMap containing the inferred "ultra paralogs" of the
     * external gene tree node with the sequence name seq_name. Sequence names
     * are the keys (String), numbers of observations are the values (Int).
     * "ultra paralogs" are to be inferred by method "inferOrthologs". Throws an
     * exception if seq_name is not found. 
     * 
     * @param seq_name
     *            sequence name of a external node of the gene trees
     * @return HashMap containing the inferred ultra paralogs
     *         (name(String)->value(Int))
     */
    public final HashMap<String, Integer> getInferredUltraParalogs( final String seq_name ) {
        if ( _up_hash_maps == null ) {
            return null;
        }
        return _up_hash_maps.get( seq_name );
    }

    /**
     * Returns the time (in ms) needed to run "inferOrthologs". Final variable
     * TIME needs to be set to true.
     * 
     * @return time (in ms) needed to run method "inferOrthologs"
     */
    public long getTime() {
        return _time;
    }

    /**
     * Infers the orthologs (as well the "super orthologs", the "subtree
     * neighbors", and the "ultra paralogs") for each external node of the gene
     * Trees in multiple tree File gene_trees_file (=output of PHYLIP NEIGHBOR,
     * for example). Tallies how many times each sequence is (super-)
     * orthologous towards the query. Tallies how many times each sequence is
     * ultra paralogous towards the query. Tallies how many times each sequence
     * is a subtree neighbor of the query. Gene duplications are inferred using
     * SDI. Modifies its argument species_tree. Is a little faster than
     * "inferOrthologs(File,Phylogeny)" since orthologs are only inferred for
     * query.
     * <p>
     * To obtain the results use the methods listed below.
     * 
     * @param gene_trees_file
     *            a File containing gene Trees in NH format, which is the result
     *            of performing a bootstrap analysis in PHYLIP
     * @param species_tree
     *            a species Phylogeny, which has species names in its species
     *            fields
     * @param query
     *            the sequence name of the squence whose orthologs are to be
     *            inferred
     */
    public void inferOrthologs( final File gene_trees_file, final Phylogeny species_tree, final String query )
            throws IOException {
        int bs = 0;
        if ( RIO.TIME ) {
            _time = System.currentTimeMillis();
        }
        if ( !gene_trees_file.exists() ) {
            throw new IllegalArgumentException( gene_trees_file.getAbsolutePath() + " does not exist." );
        }
        else if ( !gene_trees_file.isFile() ) {
            throw new IllegalArgumentException( gene_trees_file.getAbsolutePath() + " is not a file." );
        }
        // Read in first tree to get its sequence names
        // and strip species_tree.
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final Phylogeny gene_tree = factory.create( gene_trees_file, new PhyloXmlParser() )[ 0 ];
        // Removes from species_tree all species not found in gene_tree.
        PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( gene_tree, species_tree );
        // Removes from gene_tree all species not found in species_tree.
        PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gene_tree );
        _seq_names = getAllExternalSequenceNames( gene_tree );
        if ( ( _seq_names == null ) || ( _seq_names.length < 1 ) ) {
            return;
        }
        _o_hash_maps = new HashMap<String, HashMap<String, Integer>>();
        _so_hash_maps = new HashMap<String, HashMap<String, Integer>>();
        _up_hash_maps = new HashMap<String, HashMap<String, Integer>>();
        _sn_hash_maps = new HashMap<String, HashMap<String, Integer>>();
        _o_hash_maps.put( query, new HashMap<String, Integer>( _seq_names.length ) );
        _so_hash_maps.put( query, new HashMap<String, Integer>( _seq_names.length ) );
        _up_hash_maps.put( query, new HashMap<String, Integer>( _seq_names.length ) );
        _sn_hash_maps.put( query, new HashMap<String, Integer>( _seq_names.length ) );
        // Go through all gene trees in the file.
        final Phylogeny[] gene_trees = factory.create( gene_trees_file, new PhyloXmlParser() );
        for( final Phylogeny gt : gene_trees ) {
            bs++;
            // Removes from gene_tree all species not found in species_tree.
            PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gt );
            inferOrthologsHelper( gt, species_tree, query );
            // System.out.println( bs );
        }
        setBootstraps( bs );
        if ( RIO.TIME ) {
            _time = ( System.currentTimeMillis() - _time );
        }
    }

    // Helper method which performs the actual ortholog inference for
    // the external node with seqname query.
    private void inferOrthologsHelper( final Phylogeny gene_tree, final Phylogeny species_tree, final String query ) {
        Phylogeny assigned_tree = null;
        List<PhylogenyNode> nodes = null;
        final SDIR sdiunrooted = new SDIR();
        List<PhylogenyNode> orthologs = null;
        List<PhylogenyNode> super_orthologs = null;
        List<PhylogenyNode> ultra_paralogs = null;
        List<PhylogenyNode> subtree_neighbors = null;
        assigned_tree = sdiunrooted.infer( gene_tree,
                                           species_tree,
                                           RIO.ROOT_BY_MINIMIZING_MAPPING_COST,
                                           RIO.ROOT_BY_MINIMIZING_SUM_OF_DUPS,
                                           RIO.ROOT_BY_MINIMIZING_TREE_HEIGHT,
                                           true,
                                           1 )[ 0 ];
        setExtNodesOfAnalyzedGeneTrees( assigned_tree.getNumberOfExternalNodes() );
        nodes = assigned_tree.getNodesViaSequenceName( query );
        if ( nodes.size() > 1 ) {
            throw new IllegalArgumentException( "node named [" + query + "] not unique" );
        }
        else if ( nodes.isEmpty() ) {
            throw new IllegalArgumentException( "no node containing a sequence named [" + query + "] found" );
        }
        final PhylogenyNode query_node = nodes.get( 0 );
        final PhylogenyMethods methods = PhylogenyMethods.getInstance();
        orthologs = methods.getOrthologousNodes( assigned_tree, query_node );
        updateHash( _o_hash_maps, query, orthologs );
        super_orthologs = PhylogenyMethods.getSuperOrthologousNodes( query_node );
        updateHash( _so_hash_maps, query, super_orthologs );
        subtree_neighbors = getSubtreeNeighbors( query_node, 2 );
        updateHash( _sn_hash_maps, query, subtree_neighbors );
        ultra_paralogs = PhylogenyMethods.getUltraParalogousNodes( query_node );
        updateHash( _up_hash_maps, query, ultra_paralogs );
    }

    /**
     * Returns an ArrayList containg the names of orthologs of the PhylogenyNode
     * with seq name seq_name.
     * 
     * @param seq_name
     *            sequence name of a external node of the gene trees
     * @param threshold_orthologs
     *            the minimal number of observations for a a sequence to be
     *            reported as orthologous as percentage (0.0-100.0%)
     * @return ArrayList containg the names of orthologs of the PhylogenyNode
     *         with seq name seq_name
     */
    public ArrayList<String> inferredOrthologsToArrayList( final String seq_name, double threshold_orthologs ) {
        HashMap<String, Integer> o_hashmap = null;
        String name = null;
        double o = 0.0;
        final ArrayList<String> arraylist = new ArrayList<String>();
        if ( _o_hash_maps == null ) {
            throw new IllegalStateException( "Orthologs have not been calculated (successfully)." );
        }
        if ( threshold_orthologs < 0.0 ) {
            threshold_orthologs = 0.0;
        }
        else if ( threshold_orthologs > 100.0 ) {
            threshold_orthologs = 100.0;
        }
        o_hashmap = getInferredOrthologs( seq_name );
        if ( o_hashmap == null ) {
            throw new IllegalStateException( "Orthologs for " + seq_name + " were not established." );
        }
        if ( _seq_names.length > 0 ) {
            I: for( int i = 0; i < _seq_names.length; ++i ) {
                name = _seq_names[ i ];
                if ( name.equals( seq_name ) ) {
                    continue I;
                }
                o = getBootstrapValueFromHash( o_hashmap, name );
                if ( o < threshold_orthologs ) {
                    continue I;
                }
                arraylist.add( name );
            }
        }
        return arraylist;
    }

    /**
     * Returns a String containg the names of orthologs of the PhylogenyNode
     * with seq name query_name. The String also contains how many times a
     * particular ortholog has been observed.
     * <p>
     * <ul>
     * The output order is (per line): Name, Ortholog, Subtree neighbor, Super
     * ortholog, Distance
     * </ul>
     * <p>
     * The sort priority of this is determined by sort in the following manner:
     * <ul>
     * <li>0 : Ortholog
     * <li>1 : Ortholog, Super ortholog
     * <li>2 : Super ortholog, Ortholog
     * <li>3 : Ortholog, Distance
     * <li>4 : Distance, Ortholog
     * <li>5 : Ortholog, Super ortholog, Distance
     * <li>6 : Ortholog, Distance, Super ortholog
     * <li>7 : Super ortholog, Ortholog, Distance
     * <li>8 : Super ortholog, Distance, Ortholog
     * <li>9 : Distance, Ortholog, Super ortholog
     * <li>10 : Distance, Super ortholog, Ortholog
     * <li>11 : Ortholog, Subtree neighbor, Distance
     * <li>12 : Ortholog, Subtree neighbor, Super ortholog, Distance (default)
     * <li>13 : Ortholog, Super ortholog, Subtree neighbor, Distance
     * <li>14 : Subtree neighbor, Ortholog, Super ortholog, Distance
     * <li>15 : Subtree neighbor, Distance, Ortholog, Super ortholog
     * <li>16 : Ortholog, Distance, Subtree neighbor, Super ortholog
     * <li>17 : Ortholog, Subtree neighbor, Distance, Super ortholog
     * </ul>
     * <p>
     * Returns "-" if no putative orthologs have been found (given
     * threshold_orthologs).
     * <p>
     * Orthologs are to be inferred by method "inferOrthologs".
     * <p>
     * (Last modified: 05/08/01)
     * 
     * @param query_name
     *            sequence name of a external node of the gene trees
     * @param sort
     *            order and sort priority
     * @param threshold_orthologs
     *            the minimal number of observations for a a sequence to be
     *            reported as orthologous, in percents (0.0-100.0%)
     * @param threshold_subtreeneighborings
     *            the minimal number of observations for a a sequence to be
     *            reported as orthologous, in percents (0.0-100.0%)
     * @return String containing the inferred orthologs, String containing "-"
     *         if no orthologs have been found null in case of error
     * @see #inferOrthologs(File,Phylogeny,String)
     * @see #inferOrthologs(Phylogeny[],Phylogeny)
     * @see #inferOrthologs(File,Phylogeny)
     * @see #getOrder(int)
     */
    public StringBuffer inferredOrthologsToString( final String query_name,
                                                   int sort,
                                                   double threshold_orthologs,
                                                   double threshold_subtreeneighborings ) {
        HashMap<String, Integer> o_hashmap = null;
        HashMap<String, Integer> s_hashmap = null;
        HashMap<String, Integer> n_hashmap = null;
        String name = "";
        double o = 0.0, // Orthologs.
        s = 0.0, // Super orthologs.
        sn = 0.0, // Subtree neighbors.
        value1 = 0.0, value2 = 0.0, value3 = 0.0, value4 = 0.0, d = 0.0;
        final ArrayList<Tuplet> nv = new ArrayList<Tuplet>();
        if ( ( _o_hash_maps == null ) || ( _so_hash_maps == null ) || ( _sn_hash_maps == null ) ) {
            throw new IllegalStateException( "Orthologs have not been calculated (successfully)" );
        }
        if ( ( sort < 0 ) || ( sort > 17 ) ) {
            sort = 12;
        }
        if ( ( sort > 2 ) && ( _m == null ) && ( _l == null ) ) {
            throw new IllegalStateException( "Distance list or matrix have not been read in (successfully)" );
        }
        if ( threshold_orthologs < 0.0 ) {
            threshold_orthologs = 0.0;
        }
        else if ( threshold_orthologs > 100.0 ) {
            threshold_orthologs = 100.0;
        }
        if ( threshold_subtreeneighborings < 0.0 ) {
            threshold_subtreeneighborings = 0.0;
        }
        else if ( threshold_subtreeneighborings > 100.0 ) {
            threshold_subtreeneighborings = 100.0;
        }
        o_hashmap = getInferredOrthologs( query_name );
        s_hashmap = getInferredSuperOrthologs( query_name );
        n_hashmap = getInferredSubtreeNeighbors( query_name );
        if ( ( o_hashmap == null ) || ( s_hashmap == null ) || ( n_hashmap == null ) ) {
            throw new IllegalStateException( "Orthologs for " + query_name + " were not established" );
        }
        final StringBuffer orthologs = new StringBuffer();
        if ( _seq_names.length > 0 ) {
            I: for( int i = 0; i < _seq_names.length; ++i ) {
                name = _seq_names[ i ];
                if ( name.equals( query_name ) ) {
                    continue I;
                }
                o = getBootstrapValueFromHash( o_hashmap, name );
                if ( o < threshold_orthologs ) {
                    continue I;
                }
                sn = getBootstrapValueFromHash( n_hashmap, name );
                if ( sn < threshold_subtreeneighborings ) {
                    continue I;
                }
                s = getBootstrapValueFromHash( s_hashmap, name );
                if ( sort >= 3 ) {
                    if ( _m != null ) {
                        d = getDistance( query_name, name );
                    }
                    else {
                        d = getDistance( name );
                    }
                }
                switch ( sort ) {
                    case 0:
                        nv.add( new Tuplet( name, o, 5 ) );
                        break;
                    case 1:
                        nv.add( new Tuplet( name, o, s, 5 ) );
                        break;
                    case 2:
                        nv.add( new Tuplet( name, s, o, 5 ) );
                        break;
                    case 3:
                        nv.add( new Tuplet( name, o, d, 1 ) );
                        break;
                    case 4:
                        nv.add( new Tuplet( name, d, o, 0 ) );
                        break;
                    case 5:
                        nv.add( new Tuplet( name, o, s, d, 2 ) );
                        break;
                    case 6:
                        nv.add( new Tuplet( name, o, d, s, 1 ) );
                        break;
                    case 7:
                        nv.add( new Tuplet( name, s, o, d, 2 ) );
                        break;
                    case 8:
                        nv.add( new Tuplet( name, s, d, o, 1 ) );
                        break;
                    case 9:
                        nv.add( new Tuplet( name, d, o, s, 0 ) );
                        break;
                    case 10:
                        nv.add( new Tuplet( name, d, s, o, 0 ) );
                        break;
                    case 11:
                        nv.add( new Tuplet( name, o, sn, d, 2 ) );
                        break;
                    case 12:
                        nv.add( new Tuplet( name, o, sn, s, d, 3 ) );
                        break;
                    case 13:
                        nv.add( new Tuplet( name, o, s, sn, d, 3 ) );
                        break;
                    case 14:
                        nv.add( new Tuplet( name, sn, o, s, d, 3 ) );
                        break;
                    case 15:
                        nv.add( new Tuplet( name, sn, d, o, s, 1 ) );
                        break;
                    case 16:
                        nv.add( new Tuplet( name, o, d, sn, s, 1 ) );
                        break;
                    case 17:
                        nv.add( new Tuplet( name, o, sn, d, s, 2 ) );
                        break;
                    default:
                        nv.add( new Tuplet( name, o, 5 ) );
                }
            } // End of I for loop.
            if ( ( nv != null ) && ( nv.size() > 0 ) ) {
                orthologs.append( "[seq name]\t\t[ortho]\t[st-n]\t[sup-o]\t[dist]" + ForesterUtil.LINE_SEPARATOR );
                final Tuplet[] nv_array = new Tuplet[ nv.size() ];
                for( int j = 0; j < nv.size(); ++j ) {
                    nv_array[ j ] = nv.get( j );
                }
                Arrays.sort( nv_array );
                for( int i = 0; i < nv_array.length; ++i ) {
                    name = nv_array[ i ].getKey();
                    value1 = nv_array[ i ].getValue1();
                    value2 = nv_array[ i ].getValue2();
                    value3 = nv_array[ i ].getValue3();
                    value4 = nv_array[ i ].getValue4();
                    orthologs.append( addNameAndValues( name, value1, value2, value3, value4, sort ) );
                }
            }
        }
        // No orthologs found.
        if ( ( orthologs == null ) || ( orthologs.length() < 1 ) ) {
            orthologs.append( "-" );
        }
        return orthologs;
    } // inferredOrthologsToString( String, int, double )

    // Helper method for inferredOrthologTableToFile.
    // Returns individual rows for the table as String.
    private String inferredOrthologsToTableHelper( final String name2,
                                                   final String[] names,
                                                   final int j,
                                                   final boolean super_orthologs ) {
        HashMap<String, Integer> hashmap = null;
        String name = null, orthologs = new String( "" );
        int value = 0;
        if ( !super_orthologs ) {
            hashmap = getInferredOrthologs( name2 );
        }
        else {
            hashmap = getInferredSuperOrthologs( name2 );
        }
        if ( hashmap == null ) {
            throw new RuntimeException( "Unexpected failure in method inferredOrthologsToTableHelper" );
        }
        for( int i = 0; i < names.length; ++i ) {
            name = names[ i ];
            if ( !hashmap.containsKey( name ) ) {
                value = 0;
            }
            else {
                value = hashmap.get( name );
            }
            if ( i == j ) {
                // Sanity check.
                if ( value != 0 ) {
                    throw new RuntimeException( "Failed sanity check in method inferredOrthologsToTableHelper: value not 0." );
                }
                orthologs += ( " " + "\t" );
            }
            else {
                orthologs += ( value + "\t" );
            }
        }
        return orthologs;
    }

    /**
     * Writes the orthologs for each external node of the gene trees to outfile
     * in the form of a table. Orthologs are to be inferred by method
     * "inferOrthologs". Overwrites without asking! (Last modified: 12/07/00)
     * 
     * @param outfile
     *            the File to write to
     */
    public void inferredOrthologTableToFile( final File outfile ) throws IOException {
        if ( _o_hash_maps == null ) {
            return;
        }
        inferredOrthologTableToFile( outfile, false );
    }

    // Helper for inferredOrthologTableToFile(File).
    // (Last modified: 11/28/00)
    private void inferredOrthologTableToFile( final File outfile, final boolean super_orthologs ) throws IOException {
        String name = "", line = "";
        PrintWriter out = null;
        if ( _seq_names == null ) {
            throw new IllegalStateException( "inferredOrthologTableToFile: seq_names_ is null." );
        }
        Arrays.sort( _seq_names );
        out = new PrintWriter( new FileWriter( outfile ), true );
        if ( out == null ) {
            throw new RuntimeException( "inferredOrthologTableToFile: failure to create PrintWriter." );
        }
        line = "\t\t\t\t";
        for( int i = 0; i < _seq_names.length; ++i ) {
            line += ( i + ")\t" );
        }
        line += "\n";
        out.println( line );
        for( int i = 0; i < _seq_names.length; ++i ) {
            name = _seq_names[ i ];
            if ( name.length() < 8 ) {
                line = i + ")\t" + name + "\t\t\t";
            }
            else if ( name.length() < 16 ) {
                line = i + ")\t" + name + "\t\t";
            }
            else {
                line = i + ")\t" + name + "\t";
            }
            line += inferredOrthologsToTableHelper( name, _seq_names, i, super_orthologs );
            out.println( line );
        }
        out.close();
    }

    /**
     * Writes the "super orthologs" for each external nodes of the gene trees to
     * outfile in the form of a table. Super orthologs are to be inferred by
     * method "inferOrthologs". Overwrites without asking!
     * 
     * @param outfile
     *            the File to write to
     */
    public void inferredSuperOrthologTableToFile( final File outfile ) throws IOException {
        if ( _so_hash_maps == null ) {
            return;
        }
        inferredOrthologTableToFile( outfile, true );
    }

    /**
     * Returns a String containg the names of orthologs of the PhylogenyNode
     * with seq name query_name. The String also contains how many times a
     * particular ortholog has been observed. Returns "-" if no putative
     * orthologs have been found (given threshold_orthologs).
     * <p>
     * Orthologs are to be inferred by method "inferOrthologs".
     * 
     * @param query_name
     *            sequence name of a external node of the gene trees
     * @param return_dists
     * @param threshold_ultra_paralogs
     *            between 1 and 100
     * @return String containing the inferred orthologs, String containing "-"
     *         if no orthologs have been found null in case of error
     */
    public String inferredUltraParalogsToString( final String query_name,
                                                 final boolean return_dists,
                                                 double threshold_ultra_paralogs ) {
        HashMap<String, Integer> sp_hashmap = null;
        String name = "", ultra_paralogs = "";
        int sort = 0;
        double sp = 0.0, value1 = 0.0, value2 = 0.0, d = 0.0;
        final List<Tuplet> nv = new ArrayList<Tuplet>();
        if ( threshold_ultra_paralogs < 1.0 ) {
            threshold_ultra_paralogs = 1.0;
        }
        else if ( threshold_ultra_paralogs > 100.0 ) {
            threshold_ultra_paralogs = 100.0;
        }
        if ( _up_hash_maps == null ) {
            throw new IllegalStateException( "Ultra paralogs have not been calculated (successfully)." );
        }
        if ( return_dists && ( _m == null ) && ( _l == null ) ) {
            throw new IllegalStateException( "Distance list or matrix have not been read in (successfully)." );
        }
        sp_hashmap = getInferredUltraParalogs( query_name );
        if ( sp_hashmap == null ) {
            throw new IllegalStateException( "Ultra paralogs for " + query_name + " were not established" );
        }
        if ( _seq_names.length > 0 ) {
            I: for( int i = 0; i < _seq_names.length; ++i ) {
                name = _seq_names[ i ];
                if ( name.equals( query_name ) ) {
                    continue I;
                }
                sp = getBootstrapValueFromHash( sp_hashmap, name );
                if ( sp < threshold_ultra_paralogs ) {
                    continue I;
                }
                if ( return_dists ) {
                    if ( _m != null ) {
                        d = getDistance( query_name, name );
                    }
                    else {
                        d = getDistance( name );
                    }
                    nv.add( new Tuplet( name, sp, d, 1 ) );
                }
                else {
                    nv.add( new Tuplet( name, sp, 5 ) );
                }
            } // End of I for loop.
            if ( ( nv != null ) && ( nv.size() > 0 ) ) {
                final Tuplet[] nv_array = new Tuplet[ nv.size() ];
                for( int j = 0; j < nv.size(); ++j ) {
                    nv_array[ j ] = nv.get( j );
                }
                Arrays.sort( nv_array );
                if ( return_dists ) {
                    sort = 91;
                }
                else {
                    sort = 90;
                }
                for( int i = 0; i < nv_array.length; ++i ) {
                    name = nv_array[ i ].getKey();
                    value1 = nv_array[ i ].getValue1();
                    value2 = nv_array[ i ].getValue2();
                    ultra_paralogs += addNameAndValues( name, value1, value2, 0.0, 0.0, sort );
                }
            }
        }
        // No ultra paralogs found.
        if ( ( ultra_paralogs == null ) || ( ultra_paralogs.length() < 1 ) ) {
            ultra_paralogs = "-";
        }
        return ultra_paralogs;
    }

    public final void readDistanceMatrix( final File matrix_file ) throws IOException {
        DistanceMatrix[] matrices = null;
        final SymmetricalDistanceMatrixParser parser = SymmetricalDistanceMatrixParser.createInstance();
        matrices = parser.parse( matrix_file );
        if ( ( matrices == null ) || ( matrices.length == 0 ) ) {
            throw new IOException( "failed to parse distance matrix from [" + matrix_file + "]" );
        }
        if ( matrices.length > 1 ) {
            throw new IOException( "[" + matrix_file + "] contains more than once distance matrix" );
        }
        _m = matrices[ 0 ];
    }

    /**
     * Brings this into the same state as immediately after construction.
     */
    private final void reset() {
        _o_hash_maps = null;
        _so_hash_maps = null;
        _up_hash_maps = null;
        _seq_names = null;
        _m = null;
        _l = null;
        _bootstraps = 1;
        _ext_nodes_ = 0;
        _time = 0;
    }

    /**
     * Sets the numbers of trees analyzed.
     * @param the
     *            numbers of trees analyzed
     */
    private void setBootstraps( int i ) {
        if ( i < 1 ) {
            i = 1;
        }
        _bootstraps = i;
    }

    /**
     * Sets number of ext nodes in gene trees analyzed (after stripping).
     * @param the
     *            number of ext nodes in gene trees analyzed (after stripping)
     */
    private void setExtNodesOfAnalyzedGeneTrees( int i ) {
        if ( i < 1 ) {
            i = 0;
        }
        _ext_nodes_ = i;
    }

    // Helper for doInferOrthologs( Phylogeny, Phylogeny, String )
    // and doInferOrthologs( Phylogeny, Phylogeny ).
    private void updateHash( final HashMap<String, HashMap<String, Integer>> counter_map,
                             final String query_seq_name,
                             final List<PhylogenyNode> nodes ) {
        final HashMap<String, Integer> hash_map = counter_map.get( query_seq_name );
        if ( hash_map == null ) {
            throw new RuntimeException( "Unexpected failure in method updateHash." );
        }
        for( int j = 0; j < nodes.size(); ++j ) {
            final String seq_name = ( nodes.get( j ) ).getNodeData().getSequence().getName();
            if ( hash_map.containsKey( seq_name ) ) {
                hash_map.put( seq_name, hash_map.get( seq_name ) + 1 );
            }
            else {
                hash_map.put( seq_name, 1 );
            }
        }
    }

    // Helper method for inferredOrthologsToString
    // and inferredUltraParalogsToString.
    private final static String addNameAndValues( final String name,
                                                  final double value1,
                                                  final double value2,
                                                  final double value3,
                                                  final double value4,
                                                  final int sort ) {
        final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.#####" );
        df.setDecimalSeparatorAlwaysShown( false );
        String line = "";
        if ( name.length() < 8 ) {
            line += ( name + "\t\t\t" );
        }
        else if ( name.length() < 16 ) {
            line += ( name + "\t\t" );
        }
        else {
            line += ( name + "\t" );
        }
        switch ( sort ) {
            case 0:
                line += addToLine( value1, df );
                line += "-\t";
                line += "-\t";
                line += "-\t";
                break;
            case 1:
                line += addToLine( value1, df );
                line += "-\t";
                line += addToLine( value2, df );
                line += "-\t";
                break;
            case 2:
                line += addToLine( value2, df );
                line += "-\t";
                line += addToLine( value1, df );
                line += "-\t";
                break;
            case 3:
                line += addToLine( value1, df );
                line += "-\t";
                line += "-\t";
                line += addToLine( value2, df );
                break;
            case 4:
                line += addToLine( value2, df );
                line += "-\t";
                line += "-\t";
                line += addToLine( value1, df );
                break;
            case 5:
                line += addToLine( value1, df );
                line += "-\t";
                line += addToLine( value2, df );
                line += addToLine( value3, df );
                break;
            case 6:
                line += addToLine( value1, df );
                line += "-\t";
                line += addToLine( value3, df );
                line += addToLine( value2, df );
                break;
            case 7:
                line += addToLine( value2, df );
                line += "-\t";
                line += addToLine( value1, df );
                line += addToLine( value3, df );
                break;
            case 8:
                line += addToLine( value3, df );
                line += "-\t";
                line += addToLine( value1, df );
                line += addToLine( value2, df );
                break;
            case 9:
                line += addToLine( value2, df );
                line += "-\t";
                line += addToLine( value3, df );
                line += addToLine( value1, df );
                break;
            case 10:
                line += addToLine( value3, df );
                line += "-\t";
                line += addToLine( value2, df );
                line += addToLine( value1, df );
                break;
            case 11:
                line += addToLine( value1, df );
                line += addToLine( value2, df );
                line += "-\t";
                line += addToLine( value3, df );
                break;
            case 12:
                line += addToLine( value1, df );
                line += addToLine( value2, df );
                line += addToLine( value3, df );
                line += addToLine( value4, df );
                break;
            case 13:
                line += addToLine( value1, df );
                line += addToLine( value3, df );
                line += addToLine( value2, df );
                line += addToLine( value4, df );
                break;
            case 14:
                line += addToLine( value2, df );
                line += addToLine( value1, df );
                line += addToLine( value3, df );
                line += addToLine( value4, df );
                break;
            case 15:
                line += addToLine( value3, df );
                line += addToLine( value1, df );
                line += addToLine( value4, df );
                line += addToLine( value2, df );
                break;
            case 16:
                line += addToLine( value1, df );
                line += addToLine( value3, df );
                line += addToLine( value4, df );
                line += addToLine( value2, df );
                break;
            case 17:
                line += addToLine( value1, df );
                line += addToLine( value2, df );
                line += addToLine( value4, df );
                line += addToLine( value3, df );
                break;
            case 90:
                line += addToLine( value1, df );
                line += "-\t";
                break;
            case 91:
                line += addToLine( value1, df );
                line += addToLine( value2, df );
                break;
        }
        line += ForesterUtil.LINE_SEPARATOR;
        return line;
    }

    // Helper for addNameAndValues.
    private final static String addToLine( final double value, final java.text.DecimalFormat df ) {
        String s = "";
        if ( value != Tuplet.DEFAULT ) {
            s = df.format( value ) + "\t";
        }
        else {
            s = "-\t";
        }
        return s;
    }

    private static String[] getAllExternalSequenceNames( final Phylogeny phy ) {
        if ( phy.isEmpty() ) {
            return null;
        }
        int i = 0;
        final String[] names = new String[ phy.getNumberOfExternalNodes() ];
        for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
            names[ i++ ] = iter.next().getNodeData().getSequence().getName();
        }
        return names;
    }

    /**
     * Returns the order in which ortholog (o), "super ortholog" (s) and
     * distance (d) are returned and sorted (priority of sort always goes from
     * left to right), given sort. For the meaning of sort
     * 
     * @see #inferredOrthologsToString(String,int,double,double)
     *      
     * @param sort
     *            determines order and sort priority
     * @return String indicating the order
     */
    public final static String getOrder( final int sort ) {
        String order = "";
        switch ( sort ) {
            case 0:
                order = "orthologies";
                break;
            case 1:
                order = "orthologies > super orthologies";
                break;
            case 2:
                order = "super orthologies > orthologies";
                break;
            case 3:
                order = "orthologies > distance to query";
                break;
            case 4:
                order = "distance to query > orthologies";
                break;
            case 5:
                order = "orthologies > super orthologies > distance to query";
                break;
            case 6:
                order = "orthologies > distance to query > super orthologies";
                break;
            case 7:
                order = "super orthologies > orthologies > distance to query";
                break;
            case 8:
                order = "super orthologies > distance to query > orthologies";
                break;
            case 9:
                order = "distance to query > orthologies > super orthologies";
                break;
            case 10:
                order = "distance to query > super orthologies > orthologies";
                break;
            case 11:
                order = "orthologies > subtree neighbors > distance to query";
                break;
            case 12:
                order = "orthologies > subtree neighbors > super orthologies > distance to query";
                break;
            case 13:
                order = "orthologies > super orthologies > subtree neighbors > distance to query";
                break;
            case 14:
                order = "subtree neighbors > orthologies > super orthologies > distance to query";
                break;
            case 15:
                order = "subtree neighbors > distance to query > orthologies > super orthologies";
                break;
            case 16:
                order = "orthologies > distance to query > subtree neighbors > super orthologies";
                break;
            case 17:
                order = "orthologies > subtree neighbors > distance to query > super orthologies";
                break;
            default:
                order = "orthologies";
                break;
        }
        return order;
    }

    public final static StringBuffer getOrderHelp() {
        final StringBuffer sb = new StringBuffer();
        sb.append( "  0: orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  1: orthologies > super orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  2: super orthologies > orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  3: orthologies > distance to query" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  4: distance to query > orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  5: orthologies > super orthologies > distance to query" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  6: orthologies > distance to query > super orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  7: super orthologies > orthologies > distance to query" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  8: super orthologies > distance to query > orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  9: distance to query > orthologies > super orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 10: distance to query > super orthologies > orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 11: orthologies > subtree neighbors > distance to query" + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 12: orthologies > subtree neighbors > super orthologies > distance to query"
                + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 13: orthologies > super orthologies > subtree neighbors > distance to query"
                + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 14: subtree neighbors > orthologies > super orthologies > distance to query"
                + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 15: subtree neighbors > distance to query > orthologies > super orthologies"
                + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 16: orthologies > distance to query > subtree neighbors > super orthologies"
                + ForesterUtil.LINE_SEPARATOR );
        sb.append( " 17: orthologies > subtree neighbors > distance to query > super orthologies"
                + ForesterUtil.LINE_SEPARATOR );
        return sb;
    }

    private final static List<PhylogenyNode> getSubtreeNeighbors( final PhylogenyNode query, final int level ) {
        PhylogenyNode node = query;
        if ( !node.isExternal() ) {
            return null;
        }
        if ( !node.isRoot() ) {
            node = node.getParent();
        }
        if ( level == 2 ) {
            if ( !node.isRoot() ) {
                node = node.getParent();
            }
        }
        else {
            throw new IllegalArgumentException( "currently only supporting level 2 subtree neighbors " );
        }
        final List<PhylogenyNode> sn = node.getAllExternalDescendants();
        sn.remove( query );
        return sn;
    }
}
