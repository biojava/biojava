// $Id: Archaeopteryx.java,v 1.54 2009/11/24 01:36:13 cmzmasek Exp $
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

import java.io.File;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.util.ForesterUtil;

// 
// java -javaagent:shiftone-jrat.jar -cp
// $HOME/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester-atv/java/forester.jar:.
// org.forester.archaeopteryx.Archaeopteryx
// -c $HOME/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester-atv/_aptx_configuration_file
//
public final class Archaeopteryx {

    private final static boolean TEST = false; //TODO remove me!

    public static MainFrame createApplication( final Phylogeny phylogeny ) {
        final Phylogeny[] phylogenies = new Phylogeny[ 1 ];
        phylogenies[ 0 ] = phylogeny;
        return createApplication( phylogenies, "", "" );
    }

    public static MainFrame createApplication( final Phylogeny[] phylogenies ) {
        return createApplication( phylogenies, "", "" );
    }

    public static MainFrame createApplication( final Phylogeny[] phylogenies,
                                               final String config_file_name,
                                               final String title ) {
        return MainFrameApplication.createInstance( phylogenies, config_file_name, title );
    }

    public static void main( final String args[] ) {
        Phylogeny[] phylogenies = null;
        String config_filename = null;
        Configuration conf = null;
        File f = null;
        try {
            int filename_index = 0;
            if ( args.length == 0 ) {
                conf = new Configuration( null, false, false );
            }
            else if ( args.length > 0 ) {
                // check for a config file
                if ( args[ 0 ].startsWith( "-c" ) ) {
                    config_filename = args[ 1 ];
                    filename_index += 2;
                }
                conf = new Configuration( config_filename, false, false );
                if ( args.length > filename_index ) {
                    f = new File( args[ filename_index ] );
                    final String err = ForesterUtil.isReadableFile( f );
                    if ( !ForesterUtil.isEmpty( err ) ) {
                        ForesterUtil.fatalError( Constants.PRG_NAME, err );
                    }
                    boolean nhx_or_nexus = false;
                    final PhylogenyParser p = ForesterUtil.createParserDependingOnFileType( f, conf
                            .isValidatePhyloXmlAgainstSchema() );
                    if ( p instanceof NHXParser ) {
                        nhx_or_nexus = true;
                        final NHXParser nhx = ( NHXParser ) p;
                        nhx.setReplaceUnderscores( conf.isReplaceUnderscoresInNhParsing() );
                        nhx.setIgnoreQuotes( conf.isNhParsingIgnoreQuotes() );
                        TAXONOMY_EXTRACTION te = TAXONOMY_EXTRACTION.NO;
                        if ( conf.isExtractPfamTaxonomyCodesInNhParsing() ) {
                            te = TAXONOMY_EXTRACTION.PFAM_STYLE_ONLY;
                        }
                        nhx.setTaxonomyExtraction( te );
                    }
                    else if ( p instanceof NexusPhylogeniesParser ) {
                        nhx_or_nexus = true;
                        final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) p;
                        nex.setReplaceUnderscores( conf.isReplaceUnderscoresInNhParsing() );
                        nex.setIgnoreQuotes( conf.isNhParsingIgnoreQuotes() );
                    }
                    else if ( p instanceof PhyloXmlParser ) {
                        MainFrameApplication.warnIfNotPhyloXmlValidation( conf );
                    }
                    phylogenies = Util.readPhylogenies( p, f );
                    if ( nhx_or_nexus && conf.isInternalNumberAreConfidenceForNhParsing() ) {
                        for( final Phylogeny phy : phylogenies ) {
                            ForesterUtil.transferInternalNodeNamesToConfidence( phy );
                        }
                    }
                }
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( Constants.PRG_NAME, "failed to start: " + e.getLocalizedMessage() );
        }
        String title = "";
        if ( f != null ) {
            title = f.getName();
        }
        try {
            String s = "";
            if ( TEST ) {
                s = "/home/czmasek/888.xml";
                if ( ForesterUtil.isReadableFile( s ) != null ) {
                    s = "/Users/zma/888.xml";
                    if ( ForesterUtil.isReadableFile( s ) != null ) {
                        s = "C:\\888.xml";
                        if ( ForesterUtil.isReadableFile( s ) != null ) {
                            s = "C:\\Documents and Settings\\czmasek\\";
                        }
                    }
                }
            }
            if ( !TEST ) {
                MainFrameApplication.createInstance( phylogenies, conf, title );
            }
            else {
                MainFrameApplication.createInstance( ParserBasedPhylogenyFactory.getInstance()
                        .create( s, new PhyloXmlParser() ), conf, title );
            }
        }
        // catch ( final IOException ex ) {
        //     ForesterUtil.fatalError( Constants.PRG_NAME, "failed to start: " + ex.getLocalizedMessage() );
        // }
        catch ( final Exception ex ) {
            Util.unexpectedException( ex );
        }
        catch ( final Error err ) {
            Util.unexpectedError( err );
        }
    }
}