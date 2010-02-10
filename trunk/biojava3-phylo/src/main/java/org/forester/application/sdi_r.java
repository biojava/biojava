// $Id: sdi_r.java,v 1.12 2009/10/30 03:00:51 cmzmasek Exp $
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
// WWW: www.phylosoft.org

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.SDIR;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class sdi_r {

    final static private String HELP_OPTION_1           = "help";
    final static private String HELP_OPTION_2           = "h";
    final static private String MIN_MAPPING_COST_OPTION = "ml";
    final static private String MIN_DUPS_OPTION         = "md";
    final static private String MIN_HEIGHT_OPTION       = "mh";
    final static private String PRG_NAME                = "sdi_r";
    final static private String PRG_VERSION             = "1.11";
    final static private String PRG_DATE                = "2009.06.19";
    final static private String E_MAIL                  = "czmasek@burnham.org";
    final static private String WWW                     = "www.phylosoft.org";   ;
    // How many resulting trees "main" should return/display.
    private final static int    TREES_TO_RETURN         = 5;

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length == 0 ) ) {
            printHelp();
            System.exit( 0 );
        }
        if ( ( args.length < 3 ) || ( cla.getNumberOfNames() != 2 ) ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( MIN_MAPPING_COST_OPTION );
        allowed_options.add( MIN_DUPS_OPTION );
        allowed_options.add( MIN_HEIGHT_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final File outfile = new File( "sdir_outfile.xml" );
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "outfile \"" + outfile + "\" already exists" );
        }
        final File gene_tree_file = cla.getFile( 0 );
        final File species_tree_file = cla.getFile( 1 );
        boolean minimize_cost = false;
        if ( cla.isOptionSet( MIN_MAPPING_COST_OPTION ) ) {
            minimize_cost = true;
        }
        boolean minimize_sum_of_dup = false;
        if ( cla.isOptionSet( MIN_DUPS_OPTION ) ) {
            minimize_sum_of_dup = true;
        }
        boolean minimize_height = false;
        if ( cla.isOptionSet( MIN_HEIGHT_OPTION ) ) {
            minimize_height = true;
        }
        int r = 0;
        Phylogeny[] gene_trees = null;
        Phylogeny species_tree = null;
        if ( minimize_cost && minimize_sum_of_dup ) {
            minimize_sum_of_dup = false;
        }
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        try {
            final PhylogenyParser pp = new PhyloXmlParser();
            species_tree = factory.create( species_tree_file, pp )[ 0 ];
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read species tree [" + species_tree_file + "]: "
                    + e.getLocalizedMessage() );
        }
        if ( !species_tree.isRooted() ) {
            ForesterUtil.fatalError( PRG_NAME, "species tree [" + species_tree_file + "] is not rooted" );
        }
        try {
            final PhylogenyParser pp = new PhyloXmlParser();
            gene_trees = factory.create( gene_tree_file, pp );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read gene trees [" + gene_tree_file + "]: "
                    + e.getLocalizedMessage() );
        }
        // Removes from gene_tree all species not found in species_tree.
        int gene_tree_counter = 0;
        final List<Phylogeny> all_result_trees = new ArrayList<Phylogeny>();
        for( final Phylogeny gene_tree : gene_trees ) {
            r = PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gene_tree );
            ForesterUtil.programMessage( PRG_NAME, "Removed " + r + " external nodes from gene tree" );
            final SDIR sdiunrooted = new SDIR();
            final long start_time = new Date().getTime();
            Phylogeny[] result_trees = null;
            try {
                result_trees = sdiunrooted.infer( gene_tree,
                                                  species_tree,
                                                  minimize_cost,
                                                  minimize_sum_of_dup,
                                                  minimize_height,
                                                  true,
                                                  sdi_r.TREES_TO_RETURN );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
            }
            final long time_req = new Date().getTime() - start_time;
            if ( minimize_cost ) {
                ForesterUtil.programMessage( PRG_NAME, "Rooted by minimizing mapping cost L" );
                if ( minimize_height ) {
                    ForesterUtil.programMessage( PRG_NAME,
                                                 "Selected tree(s) with minimal height out of resulting trees" );
                }
                ForesterUtil.programMessage( PRG_NAME, "Number differently rooted trees minimizing criterion  : "
                        + sdiunrooted.getCount() );
                ForesterUtil.programMessage( PRG_NAME, "Minimal cost                                          : "
                        + sdiunrooted.getMinimalMappingCost() );
                ForesterUtil.programMessage( PRG_NAME, "Minimal duplications                                  : "
                        + sdiunrooted.getMinimalDuplications() );
                if ( minimize_height ) {
                    ForesterUtil.programMessage( PRG_NAME, "Phylogeny height                                      : "
                            + ForesterUtil.FORMATTER_06.format( sdiunrooted.getMinimalTreeHeight() ) );
                    ForesterUtil.programMessage( PRG_NAME, "Difference in subtree heights                         : "
                            + ForesterUtil.FORMATTER_06.format( sdiunrooted.getMinimalDiffInSubTreeHeights() ) );
                }
            }
            else if ( minimize_sum_of_dup ) {
                ForesterUtil.programMessage( PRG_NAME, "Rooted by minimizing sum of duplications" );
                if ( minimize_height ) {
                    ForesterUtil.programMessage( PRG_NAME,
                                                 "Selected tree(s) with minimal height out of resulting trees" );
                }
                ForesterUtil.programMessage( PRG_NAME, "Number differently rooted trees minimizing criterion        : "
                        + sdiunrooted.getCount() );
                ForesterUtil.programMessage( PRG_NAME, "Minimal duplications                                        : "
                        + sdiunrooted.getMinimalDuplications() );
                if ( minimize_height ) {
                    ForesterUtil.programMessage( PRG_NAME,
                                                 "Phylogeny height                                            : "
                                                         + ForesterUtil.FORMATTER_06.format( sdiunrooted
                                                                 .getMinimalTreeHeight() ) );
                    ForesterUtil.programMessage( PRG_NAME,
                                                 "Difference in subtree heights                               : "
                                                         + ForesterUtil.FORMATTER_06.format( sdiunrooted
                                                                 .getMinimalDiffInSubTreeHeights() ) );
                }
            }
            else if ( minimize_height ) {
                ForesterUtil.programMessage( PRG_NAME, "Rooted by minimizing tree height (midpoint rooting)." );
                ForesterUtil.programMessage( PRG_NAME, "Minimal tree height                  : "
                        + ForesterUtil.FORMATTER_06.format( sdiunrooted.getMinimalTreeHeight() ) );
                ForesterUtil.programMessage( PRG_NAME, "Minimal difference in subtree heights: "
                        + ForesterUtil.FORMATTER_06.format( sdiunrooted.getMinimalDiffInSubTreeHeights() ) );
                ForesterUtil.programMessage( PRG_NAME, "Duplications in midpoint rooted tree : "
                        + sdiunrooted.getMinimalDuplications() );
            }
            else {
                ForesterUtil.programMessage( PRG_NAME, "No (re) rooting was performed." );
                ForesterUtil.programMessage( PRG_NAME, "Duplications in tree: " + sdiunrooted.getMinimalDuplications() );
            }
            ForesterUtil.programMessage( PRG_NAME, "Time requirement (minus I/O)                          : "
                    + time_req + "ms" );
            for( int i = 0; i < result_trees.length; ++i ) {
                final String name = result_trees[ i ].getName();
                if ( ForesterUtil.isEmpty( name ) ) {
                    result_trees[ i ].setName( "SDIR result [gene tree + " + gene_tree_counter + "]" + " " + i );
                }
                else {
                    result_trees[ i ].setName( name + " SDIR result [gene tree + " + gene_tree_counter + "]" + " " + i );
                }
                all_result_trees.add( result_trees[ i ] );
            }
            ++gene_tree_counter;
        } // for( final Phylogeny gene_tree : gene_trees ) 
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( outfile, all_result_trees, 0, ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failure to write output to [" + outfile + "]: "
                    + e.getLocalizedMessage() );
        }
        ForesterUtil.programMessage( PRG_NAME, "Wrote: " + outfile );
        ForesterUtil.programMessage( PRG_NAME, "OK." );
    }

    private static void printHelp() {
        System.out.println( "Usage: " + PRG_NAME
                + " <options> <gene tree(s) in phyloXML format> <species tree in phyloXML format>\"" );
        System.out.println( "\nOptions:" );
        System.out.println( " -" + MIN_MAPPING_COST_OPTION
                + " to root by minimizing the mapping cost L (and also the sum of duplications)" );
        System.out.println( " -" + MIN_DUPS_OPTION + " to root by minimizing the sum of duplications" );
        System.out.println( " -" + MIN_HEIGHT_OPTION
                + " to root by minimizing tree height (can be used together with -" + MIN_MAPPING_COST_OPTION + " or -"
                + MIN_DUPS_OPTION + ")" );
        System.out.println( "" );
    }
}
