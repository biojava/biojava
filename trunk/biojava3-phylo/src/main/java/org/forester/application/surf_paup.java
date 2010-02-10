// $Id: surf_paup.java,v 1.7 2008/08/20 22:32:56 cmzmasek Exp $
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
// WWW: www.phylosoft.org

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.nexus.NexusCharactersParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nexus.PaupLogParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogenyinference.CharacterStateMatrix;
import org.forester.phylogenyinference.CharacterStateMatrix.BinaryStates;
import org.forester.phylogenyinference.CharacterStateMatrix.Format;
import org.forester.surfacing.DomainParsimonyCalculator;
import org.forester.surfacing.SurfacingUtil;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class surf_paup {

    final static private String PRG_VERSION   = "0.90";
    final static private String PRG_DATE      = "2008.03.28";
    final static private String E_MAIL        = "czmasek@burnham.org";
    final static private String WWW           = "www.phylosoft.org/forester/applications/surfacing";
    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    private static final String PRG_NAME      = "surf_paup";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( HELP_OPTION_1 );
        allowed_options.add( HELP_OPTION_2 );
        if ( ( args.length < 2 ) ) {
            printHelp();
            System.exit( -1 );
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) ) {
            printHelp();
            System.exit( 0 );
        }
        if ( cla.getNumberOfNames() != 3 ) {
            printHelp();
            System.exit( -1 );
        }
        final File surfacing_nexus_outfile = cla.getFile( 0 );
        final File paup_log_file = cla.getFile( 1 );
        final String outfile_name = cla.getFile( 2 ).toString();
        final NexusCharactersParser nex_char_parser = new NexusCharactersParser();
        try {
            nex_char_parser.setSource( surfacing_nexus_outfile );
            nex_char_parser.parse();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "problem with parsing character labels from  ["
                    + surfacing_nexus_outfile + "]: " + e.getMessage() );
            e.printStackTrace();
        }
        final String[] labels = nex_char_parser.getCharStateLabels();
        ForesterUtil.programMessage( PRG_NAME, "read in " + labels.length + " character labels" );
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final NexusPhylogeniesParser phylogeny_parser = new NexusPhylogeniesParser();
        Phylogeny[] phylogenies = null;
        try {
            phylogenies = factory.create( surfacing_nexus_outfile, phylogeny_parser );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "problem with parsing phylogeny [" + surfacing_nexus_outfile + "]: "
                    + e.getMessage() );
            e.printStackTrace();
        }
        if ( phylogenies.length != 1 ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to parse one phylogeny from [" + surfacing_nexus_outfile
                    + "], got " + phylogenies.length + " instead" );
        }
        final Phylogeny phylogeny = phylogenies[ 0 ];
        if ( !phylogeny.isRooted() ) {
            ForesterUtil.fatalError( PRG_NAME, "phylogeny from [" + surfacing_nexus_outfile + "] is not rooted" );
        }
        ForesterUtil.postOrderRelabelInternalNodes( phylogeny, phylogeny.getNumberOfExternalNodes() + 1 );
        CharacterStateMatrix<BinaryStates> matrix = null;
        final PaupLogParser paup_log_parser = new PaupLogParser();
        try {
            paup_log_parser.setSource( paup_log_file );
            matrix = paup_log_parser.parse();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to parse matrix from  [" + paup_log_file + "]: "
                    + e.getMessage() );
        }
        ForesterUtil.programMessage( PRG_NAME, "read in character state matrix of size "
                + matrix.getNumberOfIdentifiers() + "x" + matrix.getNumberOfCharacters() );
        final DomainParsimonyCalculator domain_parsimony = DomainParsimonyCalculator.createInstance( phylogeny );
        domain_parsimony.executeOnGivenBinaryStatesMatrix( matrix, labels );
        final String sep = ForesterUtil.LINE_SEPARATOR + "###################" + ForesterUtil.LINE_SEPARATOR;
        SurfacingUtil.writeMatrixToFile( domain_parsimony.getGainLossMatrix(),
                                         outfile_name + "_paup_gl",
                                         Format.FORESTER );
        SurfacingUtil.writeMatrixToFile( domain_parsimony.getGainLossCountsMatrix(),
                                         outfile_name + "_paup_glc",
                                         Format.FORESTER );
        SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                           CharacterStateMatrix.GainLossStates.GAIN,
                                                           outfile_name + "_paup_gains",
                                                           sep,
                                                           ForesterUtil.LINE_SEPARATOR,
                                                           null );
        SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                           CharacterStateMatrix.GainLossStates.LOSS,
                                                           outfile_name + "_paup_losses",
                                                           sep,
                                                           ForesterUtil.LINE_SEPARATOR,
                                                           null );
        SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(), null, outfile_name
                + "_paup_present", sep, ForesterUtil.LINE_SEPARATOR, null );
        final String date_time = ForesterUtil.getCurrentDateTime();
        SurfacingUtil.preparePhylogeny( phylogeny, domain_parsimony, date_time, "parsimony (paup)", "paup_"
                + outfile_name, "" );
        SurfacingUtil.writePhylogenyToFile( phylogeny, outfile_name + "_paup.xml" );
        ForesterUtil.programMessage( PRG_NAME, "OK" );
    }

    private static void printHelp() {
        System.out.println();
        System.out.println( "Usage:" );
        System.out.println();
        System.out
                .println( "% java  -cp forester.jar org.forester.applications."
                        + PRG_NAME
                        + " <surfacing nexus outfile with character labels and tree> <paup log file with reconstructed states matrix> <outfile name base>" );
        System.out.println();
    }
}