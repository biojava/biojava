// $Id: pccx.java,v 1.11 2009/11/20 22:22:10 cmzmasek Exp $
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

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.pccx.BasicExternalNodeBasedCoverageExtender;
import org.forester.pccx.Coverage;
import org.forester.pccx.CoverageCalculationOptions;
import org.forester.pccx.CoverageCalculator;
import org.forester.pccx.CoverageExtender;
import org.forester.pccx.ExternalNodeBasedCoverageMethod;
import org.forester.pccx.ExternalNodeBasedCoverageMethodOptions;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

/*
 * @author Christian M. Zmasek
 */
public class pccx {

    final static private int    EXTEND_BY_DEFAULT                   = -100;
    final static private String HELP_OPTION_1                       = "help";
    final static private String HELP_OPTION_2                       = "h";
    final static private String USE_REAL_BL_OPTION                  = "d";
    final static private String USE_LOG_REAL_BL_OPTION              = "ld";
    final static private String EXTEND_BY_OPTION                    = "x";
    final static private String OUTPUT_OPTION                       = "o";
    final static private String INPUT_OPTION                        = "i";
    final static private String OUTPUT_ANNOTATED_PHYLOGENIES_OPTION = "p";
    final static private String PRG_NAME                            = "pccx";
    final static private String PRG_VERSION                         = "1.0.0";
    final static private String BRANCH_LENGTH_BASED_SCORING         = "org.forester.tools.modeling.BranchLengthBasedScoringMethod";
    final static private String BRANCH_COUNTING_BASED_SCORING       = "org.forester.tools.modeling.BranchCountingBasedScoringMethod";
    final static private String LOG_BRANCH_LENGTH_BASED_SCORING     = "org.forester.tools.modeling.LogBranchLengthBasedScoringMethod";
    final static private String PRG_DATE                            = "2008.03.04";
    final static private String WWW                                 = "www.phylosoft.org/forester/applications/pccx";
    final static private String E_MAIL                              = "czmasek@burnham.org";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( pccx.PRG_NAME, pccx.PRG_VERSION, pccx.PRG_DATE, pccx.E_MAIL, pccx.WWW );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( pccx.HELP_OPTION_1 ) || cla.isOptionSet( pccx.HELP_OPTION_2 ) ) {
            System.out.println();
            pccx.printHelp();
            System.exit( 0 );
        }
        if ( ( args.length < 2 ) ) {
            System.out.println();
            System.out.println( "Incorrect number of arguments." );
            System.out.println();
            pccx.printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        boolean use_bl = false;
        boolean use_log_bl = false;
        int extend_by = pccx.EXTEND_BY_DEFAULT;
        allowed_options.add( pccx.USE_REAL_BL_OPTION );
        allowed_options.add( pccx.USE_LOG_REAL_BL_OPTION );
        allowed_options.add( pccx.EXTEND_BY_OPTION );
        allowed_options.add( pccx.INPUT_OPTION );
        allowed_options.add( pccx.OUTPUT_OPTION );
        allowed_options.add( pccx.OUTPUT_ANNOTATED_PHYLOGENIES_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( pccx.PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        if ( cla.getNumberOfNames() < 1 ) {
            System.out.println();
            System.out.println( "No phylogenies infile indicated." );
            System.out.println();
            pccx.printHelp();
            System.exit( -1 );
        }
        final File phylogenies_infile = cla.getFile( 0 );
        final List<String> external_otu_names = new ArrayList<String>();
        if ( cla.getNumberOfNames() > 1 ) {
            for( int i = 1; i < cla.getNumberOfNames(); ++i ) {
                external_otu_names.add( cla.getName( i ) );
            }
        }
        if ( cla.isOptionSet( pccx.USE_REAL_BL_OPTION ) ) {
            use_bl = true;
        }
        if ( cla.isOptionSet( pccx.USE_LOG_REAL_BL_OPTION ) ) {
            use_log_bl = true;
        }
        if ( use_bl && use_log_bl ) {
            System.out.println();
            pccx.printHelp();
            System.exit( -1 );
        }
        if ( cla.isOptionSet( pccx.EXTEND_BY_OPTION ) ) {
            extend_by = 0;
            try {
                extend_by = cla.getOptionValueAsInt( pccx.EXTEND_BY_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( pccx.PRG_NAME, e.getMessage() );
            }
        }
        Phylogeny[] phylogenies = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ForesterUtil.createParserDependingOnFileType( phylogenies_infile, true );
            phylogenies = factory.create( phylogenies_infile, pp );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( pccx.PRG_NAME, "could not read \"" + phylogenies_infile + "\": " + e.getMessage() );
        }
        final List<Phylogeny> phylogenies_list = Arrays.asList( phylogenies );
        File outfile = null;
        PrintStream out = System.out;
        if ( cla.isOptionSet( pccx.OUTPUT_OPTION ) ) {
            try {
                outfile = new File( cla.getOptionValue( pccx.OUTPUT_OPTION ) );
                final String error = ForesterUtil.isWritableFile( outfile );
                if ( !ForesterUtil.isEmpty( error ) ) {
                    ForesterUtil.fatalError( pccx.PRG_NAME, error );
                }
                out = new PrintStream( outfile );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( pccx.PRG_NAME, e.getMessage() );
            }
        }
        File infile = null;
        BasicTable<String> intable = null;
        if ( cla.isOptionSet( pccx.INPUT_OPTION ) ) {
            try {
                infile = new File( cla.getOptionValue( pccx.INPUT_OPTION ) );
                final String error = ForesterUtil.isReadableFile( infile );
                if ( !ForesterUtil.isEmpty( error ) ) {
                    ForesterUtil.fatalError( pccx.PRG_NAME, error );
                }
                intable = BasicTableParser.parse( infile, " ", false );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( pccx.PRG_NAME, "failed to read \"" + infile + "\" [" + e.getMessage() + "]" );
            }
            try {
                for( int row = 0; row < intable.getNumberOfRows(); ++row ) {
                    System.out.println( "Adding external node: " + intable.getValueAsString( 0, row ) );
                    external_otu_names.add( intable.getValueAsString( 0, row ) );
                }
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( pccx.PRG_NAME, e.getMessage() );
            }
        }
        File annotated_phylogenies_outfile = null;
        boolean output_annoted_phylogenies = false;
        if ( cla.isOptionSet( pccx.OUTPUT_ANNOTATED_PHYLOGENIES_OPTION ) ) {
            output_annoted_phylogenies = true;
            annotated_phylogenies_outfile = new File( cla.getOptionValue( pccx.OUTPUT_ANNOTATED_PHYLOGENIES_OPTION ) );
            final String error = ForesterUtil.isWritableFile( annotated_phylogenies_outfile );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( pccx.PRG_NAME, error );
            }
        }
        try {
            final CoverageCalculationOptions options;
            if ( use_log_bl ) {
                options = new ExternalNodeBasedCoverageMethodOptions( pccx.LOG_BRANCH_LENGTH_BASED_SCORING );
            }
            else if ( use_bl ) {
                options = new ExternalNodeBasedCoverageMethodOptions( pccx.BRANCH_LENGTH_BASED_SCORING );
            }
            else {
                options = new ExternalNodeBasedCoverageMethodOptions( pccx.BRANCH_COUNTING_BASED_SCORING );
            }
            final int s = phylogenies_list.get( 0 ).getNumberOfExternalNodes() - external_otu_names.size();
            if ( extend_by > s ) {
                extend_by = s;
            }
            System.out.println();
            System.out.println( "Options: " + options.asString() );
            System.out.println();
            if ( extend_by != pccx.EXTEND_BY_DEFAULT ) {
                if ( extend_by > 0 ) {
                    System.out.println( "Printing " + extend_by + " names to extend coverage in an optimal manner:" );
                }
                else {
                    System.out.println( "Printing names to completely extend coverage in an optimal manner:" );
                }
                System.out.println();
                final CoverageCalculator cc = CoverageCalculator.getInstance( new ExternalNodeBasedCoverageMethod(),
                                                                              options );
                final CoverageExtender ce = new BasicExternalNodeBasedCoverageExtender();
                Coverage cov = cc.calculateCoverage( phylogenies_list, external_otu_names, false );
                System.out.println( " before:" );
                System.out.println( cov.asString() );
                System.out.println();
                final List<String> result = ce.find( phylogenies_list, external_otu_names, extend_by, options, out );
                final List<String> new_names = new ArrayList<String>( external_otu_names );
                for( final Object element : result ) {
                    final String n = ( String ) element;
                    new_names.add( n );
                }
                cov = cc.calculateCoverage( phylogenies_list, new_names, output_annoted_phylogenies );
                System.out.println();
                System.out.println( " after:" );
                System.out.println( cov.asString() );
            }
            else {
                final CoverageCalculator cc = CoverageCalculator.getInstance( new ExternalNodeBasedCoverageMethod(),
                                                                              options );
                final Coverage cov = cc.calculateCoverage( phylogenies_list,
                                                           external_otu_names,
                                                           output_annoted_phylogenies );
                System.out.println( cov.asString() );
            }
            System.out.println();
            if ( output_annoted_phylogenies ) {
                try {
                    final PhylogenyWriter writer = new PhylogenyWriter();
                    writer.toPhyloXML( annotated_phylogenies_outfile, phylogenies_list.get( 0 ), 1 );
                    System.out.println( "Wrote annotated phylogeny to \"" + annotated_phylogenies_outfile + "\"" );
                    System.out.println();
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( pccx.PRG_NAME, "Failed to write to \"" + annotated_phylogenies_outfile
                            + "\" [" + e.getMessage() + "]" );
                }
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( pccx.PRG_NAME, e.toString() );
        }
        System.out.println();
        System.out.println( pccx.PRG_NAME + ": successfully completed" );
        System.out.println( "If this application is useful to you, please cite:" );
        System.out.println( pccx.WWW );
        System.out.println();
        out.flush();
        out.close();
    }

    private static void printHelp() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( pccx.PRG_NAME
                + "  [options] <phylogen(y|ies) infile> [external node name 1] [name 2] ... [name n]" );
        System.out.println();
        System.out.println( " Options: " );
        System.out.println();
        System.out.println( " -d        : 1/distance based scoring method (instead of branch counting based)" );
        System.out.println( " -ld       : -ln(distance) based scoring method (instead of branch counting based)" );
        System.out.println( " -x[=<n>]  : optimally extend coverage by <n> external nodes. Use none, 0," );
        System.out.println( "             or negative value for complete coverage extension." );
        System.out.println( " -o=<file> : write output to <file>" );
        System.out.println( " -i=<file> : read (new-line separated) external node names from <file>" );
        System.out.println( " -" + pccx.OUTPUT_ANNOTATED_PHYLOGENIES_OPTION
                + "=<file> : write output as annotated phylogeny to <file> (only first" );
        System.out.println( "             phylogeny in phylogenies infile is used)" );
        System.out.println();
    }
}
