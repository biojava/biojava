// $Id: meta_ontologizer.java,v 1.12 2009/04/29 22:35:50 cmzmasek Exp $
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
import java.util.List;

import org.forester.go.PfamToGoMapping;
import org.forester.go.PfamToGoParser;
import org.forester.go.etc.MetaOntologizer;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class meta_ontologizer {

    final static private String HELP_OPTION_1      = "help";
    final static private String HELP_OPTION_2      = "h";
    final static private String P_OPTION           = "p";
    final static private String PRG_NAME           = "meta_ontologizer";
    final static private String PRG_VERSION        = "1.10";
    final static private String PRG_DATE           = "2009.04.29";
    final static private String E_MAIL             = "czmasek@burnham.org";
    final static private String WWW                = "www.phylosoft.org/forester/";
    private static final String RESULT_FILE_PREFIX = "table";

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
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( P_OPTION );
        final List<String> mandatory_options = new ArrayList<String>();
        mandatory_options.add( P_OPTION );
        if ( ( cla.getNumberOfNames() != 5 ) && ( cla.getNumberOfNames() != 6 ) ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final String missing = cla.validateMandatoryOptionsAsString( mandatory_options );
        if ( missing.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "missing option(s): " + missing );
        }
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final File obo_file = cla.getFile( 0 );
        final File pfam2go_file = cla.getFile( 1 );
        final File ontologizer_outdir = cla.getFile( 2 );
        File domain_gain_loss_file = null;
        String outfile_base = null;
        String comment = null;
        if ( cla.getNumberOfNames() == 6 ) {
            domain_gain_loss_file = cla.getFile( 3 );
            outfile_base = cla.getName( 4 );
            comment = cla.getName( 5 );
        }
        else {
            outfile_base = cla.getName( 3 );
            comment = cla.getName( 4 );
        }
        double p_adjusted_upper_limit = -1;
        try {
            p_adjusted_upper_limit = cla.getOptionValueAsDouble( P_OPTION );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        try {
            final PfamToGoParser parser = new PfamToGoParser( pfam2go_file );
            final List<PfamToGoMapping> pfam_to_go_mappings = parser.parse();
            ForesterUtil.programMessage( PRG_NAME, "parsed " + pfam_to_go_mappings.size() + " Pfam to GO mappings" );
            MetaOntologizer.reformat( ontologizer_outdir,
                                      RESULT_FILE_PREFIX,
                                      domain_gain_loss_file,
                                      outfile_base,
                                      obo_file,
                                      p_adjusted_upper_limit,
                                      comment,
                                      pfam_to_go_mappings );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            e.printStackTrace();
        }
        ForesterUtil.programMessage( PRG_NAME, "OK" );
        System.out.println();
    }

    private static void printHelp() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out
                .println( PRG_NAME
                        + " -p=<adj P value limit> <obo file> <pfam to go file> <ontologizer outdir> [domain gain loss file] <base for meta ontologizer outfile> <comment>" );
        System.out.println();
    }
}
