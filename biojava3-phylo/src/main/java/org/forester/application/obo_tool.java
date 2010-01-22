// $Id: obo_tool.java,v 1.3 2009/01/13 19:49:32 cmzmasek Exp $
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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.go.GoTerm;
import org.forester.go.OBOparser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class obo_tool {

    private static final String IDS_TO_NAMES_SUFFIX  = "_ids_to_names";
    final static private String HELP_OPTION_1        = "help";
    final static private String HELP_OPTION_2        = "h";
    final static private String GO_ID_TO_NAME_OPTION = "i";
    final static private String PRG_NAME             = "obo_tool";
    final static private String PRG_VERSION          = "1.00";
    final static private String PRG_DATE             = "2008.11.26";
    final static private String E_MAIL               = "czmasek@burnham.org";
    final static private String WWW                  = "www.phylosoft.org/forester/";

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
        if ( args.length < 3 ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( GO_ID_TO_NAME_OPTION );
        if ( cla.getNumberOfNames() != 2 ) {
            System.out.println();
            System.out.println( "[" + PRG_NAME + "] incorrect number of arguments" );
            System.out.println();
            printHelp();
            System.exit( -1 );
        }
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        boolean output_ids_to_names = false;
        if ( cla.isOptionSet( GO_ID_TO_NAME_OPTION ) ) {
            output_ids_to_names = true;
        }
        final File infile = cla.getFile( 0 );
        final File outfile = cla.getFile( 1 );
        final OBOparser parser = new OBOparser( infile, OBOparser.ReturnType.BASIC_GO_TERM );
        List<GoTerm> go_terms = null;
        try {
            go_terms = parser.parse();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.toString() );
        }
        ForesterUtil.programMessage( PRG_NAME, "successfully read in " + go_terms.size() + " GO terms from [" + infile
                + "]" );
        if ( output_ids_to_names ) {
            final File outfile_ids_to_names = new File( outfile + IDS_TO_NAMES_SUFFIX );
            final String error = ForesterUtil.isWritableFile( outfile_ids_to_names );
            if ( !ForesterUtil.isEmpty( error ) ) {
                ForesterUtil.fatalError( PRG_NAME, error );
            }
            try {
                final Writer out = new BufferedWriter( new FileWriter( outfile_ids_to_names ) );
                for( final GoTerm go_term : go_terms ) {
                    out.write( go_term.getGoId().getId() );
                    out.write( "\t" );
                    out.write( go_term.getDefinition() );
                    out.write( ForesterUtil.LINE_SEPARATOR );
                }
                out.close();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.toString() );
            }
            ForesterUtil.programMessage( PRG_NAME, "wrote: [" + outfile_ids_to_names + "]" );
        }
        ForesterUtil.programMessage( PRG_NAME, "OK" );
        System.out.println();
    }

    private static void printHelp() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " <options> <obo infile> <outfile>" );
        System.out.println();
        System.out.println( " options: " );
        System.out.println();
        System.out.println( "   -" + GO_ID_TO_NAME_OPTION + ": output GO id to name map file" );
        System.out.println();
    }
}
