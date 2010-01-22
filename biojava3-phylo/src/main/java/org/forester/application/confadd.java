// $Id: confadd.java,v 1.7 2009/12/17 02:28:13 cmzmasek Exp $
//
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
import java.util.ArrayList;
import java.util.List;

import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.tools.ConfidenceAssessor;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class confadd {

    final static private String HELP_OPTION_1    = "help";
    final static private String HELP_OPTION_2    = "h";
    final static private String FIRST_OPTION     = "f";
    final static private String LAST_OPTION      = "l";
    final static private String STRICT_OPTION    = "s";
    final static private String NORMALIZE_OPTION = "n";
    final static private String PRG_NAME         = "confadd";
    final static private String PRG_VERSION      = "1.00 beta 1";
    final static private String PRG_DATE         = "2009.12.16";
    final static private String E_MAIL           = "czmasek@burnham.org";
    final static private String WWW              = "www.phylosoft.org/forester/";

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
        if ( args.length < 4 ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        if ( cla.getNumberOfNames() != 4 ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( FIRST_OPTION );
        allowed_options.add( LAST_OPTION );
        allowed_options.add( STRICT_OPTION );
        allowed_options.add( NORMALIZE_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final String confidence_type = cla.getName( 0 );
        final File target_file = cla.getFile( 1 );
        final File evaluators_file = cla.getFile( 2 );
        final File outfile = cla.getFile( 3 );
        if ( ForesterUtil.isEmpty( confidence_type ) ) {
            ForesterUtil.fatalError( PRG_NAME, "attempt to use empty confidence type" );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        if ( !target_file.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "target [" + target_file + "] does not exist" );
        }
        if ( !evaluators_file.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "evaluators [" + evaluators_file + "] does not exist" );
        }
        boolean strict = false;
        int first = 0;
        int last = 0;
        double norm = 0;
        try {
            if ( cla.isOptionSet( STRICT_OPTION ) ) {
                if ( cla.isOptionHasAValue( STRICT_OPTION ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "no value allowed for -" + STRICT_OPTION + " allowed" );
                }
                strict = true;
            }
            if ( cla.isOptionSet( FIRST_OPTION ) ) {
                first = cla.getOptionValueAsInt( FIRST_OPTION );
            }
            if ( cla.isOptionSet( LAST_OPTION ) ) {
                last = cla.getOptionValueAsInt( LAST_OPTION );
            }
            if ( cla.isOptionSet( NORMALIZE_OPTION ) ) {
                norm = cla.getOptionValueAsDouble( NORMALIZE_OPTION );
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "error in command line: " + e.getLocalizedMessage() );
        }
        if ( ( first < 0 ) || ( last < 0 ) ) {
            ForesterUtil
                    .fatalError( PRG_NAME,
                                 "attempt to set first or last evaluator topology to use to a number less than zero" );
        }
        if ( norm < 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "illegal value for normalizer [" + norm + "]" );
        }
        Phylogeny[] targets = null;
        Phylogeny[] evaluators = null;
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        try {
            targets = factory.create( target_file, ForesterUtil.createParserDependingOnFileType( target_file, true ) );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read target phylogenies from [" + target_file + "]: "
                    + e.getLocalizedMessage() );
        }
        if ( targets.length == 1 ) {
            ForesterUtil.programMessage( PRG_NAME, "read in one target" );
        }
        else {
            ForesterUtil.programMessage( PRG_NAME, "read in a total of " + targets.length + " targets" );
        }
        try {
            evaluators = factory.create( evaluators_file, ForesterUtil
                    .createParserDependingOnFileType( evaluators_file, true ) );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read evaluator topologies from [" + evaluators_file + "]: "
                    + e.getLocalizedMessage() );
        }
        ForesterUtil.programMessage( PRG_NAME, "read in a total of " + evaluators.length + " evaluator topologies" );
        System.gc();
        if ( last == 0 ) {
            last = evaluators.length - 1;
        }
        if ( ( last >= evaluators.length ) || ( last <= first ) ) {
            ForesterUtil.fatalError( PRG_NAME, "illegal value for first or last evaluator topology to use [" + first
                    + ", " + last + "]" );
        }
        double value = 1;
        if ( norm > 0 ) {
            value = norm / ( 1 + last - first );
        }
        ForesterUtil.programMessage( PRG_NAME, "first topology to use: " + first );
        String is_last = "";
        if ( last == ( evaluators.length - 1 ) ) {
            is_last = " (corresponds to last topology in file)";
        }
        ForesterUtil.programMessage( PRG_NAME, "last topology to use : " + last + is_last );
        ForesterUtil.programMessage( PRG_NAME, "sum of topologies used as evaluators: " + ( last - first + 1 ) );
        if ( norm > 0 ) {
            ForesterUtil.programMessage( PRG_NAME, "normalizer: " + norm + " (" + ForesterUtil.round( value, 6 ) + ")" );
        }
        else {
            ForesterUtil.programMessage( PRG_NAME, "normalizer: n/a" );
        }
        ForesterUtil.programMessage( PRG_NAME, "strict: " + strict );
        for( final Phylogeny target : targets ) {
            try {
                ConfidenceAssessor.evaluate( confidence_type, evaluators, target, strict, value, first, last );
            }
            catch ( final IllegalArgumentException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
            }
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( targets, 0, outfile, ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outfile + "]: " + e.getLocalizedMessage() );
        }
        ForesterUtil.programMessage( PRG_NAME, "wrote output to: [" + outfile + "]" );
        ForesterUtil.programMessage( PRG_NAME, "OK" );
        System.out.println();
    }

    private static void printHelp() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME
                + " [options] <confidence type> <target tree file> <evaluators tree file> <outfile>" );
        System.out.println();
        System.out.println( "options:" );
        System.out.println();
        System.out.println( " -" + STRICT_OPTION
                + "    : strict [default: non-strict]: all nodes between 'target' and 'evaluators' must match" );
        System.out.println( " -" + NORMALIZE_OPTION
                + "=<d>: normalize to this value (e.g. 100 for most bootstrap analyses) [default: no normalization]" );
        System.out.println( " -" + FIRST_OPTION + "=<i>: first evaluator topology to use (0-based) [default: 0]" );
        System.out.println( " -" + LAST_OPTION
                + "=<i>: last evaluator topology to use (0-based) [default: use all until final topology]" );
        System.out.println();
    }
}
