// $Id: NexusPhylogeniesParser.java,v 1.15 2009/10/26 23:29:40 cmzmasek Exp $
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

package org.forester.io.parsers.nexus;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.archaeopteryx.Constants;
import org.forester.io.parsers.ParserUtils;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.nhx.NHXFormatException;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class NexusPhylogeniesParser implements PhylogenyParser {

    final private static String  begin_trees               = NexusConstants.BEGIN_TREES.toLowerCase();
    final private static String  taxlabels                 = NexusConstants.TAXLABELS.toLowerCase();
    final private static String  translate                 = NexusConstants.TRANSLATE.toLowerCase();
    final private static String  tree                      = NexusConstants.TREE.toLowerCase();
    final private static String  utree                     = NexusConstants.UTREE.toLowerCase();
    final private static String  end                       = NexusConstants.END.toLowerCase();
    final private static String  endblock                  = "endblock";
    final private static Pattern TREE_NAME_PATTERN         = Pattern.compile( "\\s*.?Tree\\s+(.+?)\\s*=.+",
                                                                              Pattern.CASE_INSENSITIVE );
    final private static Pattern ROOTEDNESS_PATTERN        = Pattern.compile( ".+=\\s*\\[&([R|U])\\].*" );
    private Object               _nexus_source;
    private List<Phylogeny>      _phylogenies;
    private List<String>         _taxlabels;
    private Map<String, String>  _translate_map;
    private boolean              _replace_underscores      = NHXParser.REPLACE_UNDERSCORES_DEFAULT;
    private boolean              _ignore_quotes_in_nh_data = Constants.NH_PARSING_IGNORE_QUOTES_DEFAULT;

    private void createPhylogeny( final String name,
                                  final StringBuffer nhx,
                                  final boolean rooted_info_present,
                                  final boolean is_rooted ) throws IOException {
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final NHXParser pars = new NHXParser();
        pars.setTaxonomyExtraction( TAXONOMY_EXTRACTION.NO );
        pars.setReplaceUnderscores( isReplaceUnderscores() );
        pars.setIgnoreQuotes( isIgnoreQuotes() );
        if ( rooted_info_present ) {
            pars.setGuessRootedness( false );
        }
        final Phylogeny p = factory.create( nhx, pars )[ 0 ];
        p.setName( name );
        if ( rooted_info_present ) {
            p.setRooted( is_rooted );
        }
        if ( ( getTaxlabels().size() > 0 ) || ( getTranslateMap().size() > 0 ) ) {
            final PhylogenyNodeIterator it = p.iteratorExternalForward();
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                if ( ( getTranslateMap().size() > 0 ) && getTranslateMap().containsKey( node.getNodeName() ) ) {
                    node.setName( getTranslateMap().get( node.getNodeName() ).replaceAll( "['\"]+", "" ) );
                }
                else if ( getTaxlabels().size() > 0 ) {
                    int i = -1;
                    try {
                        i = Integer.parseInt( node.getNodeName() );
                    }
                    catch ( final NumberFormatException e ) {
                        // Ignore.
                    }
                    if ( i > 0 ) {
                        node.setName( getTaxlabels().get( i - 1 ).replaceAll( "['\"]+", "" ) );
                    }
                }
            }
        }
        getPhylogenies().add( p );
    }

    private Object getNexusSource() {
        return _nexus_source;
    }

    private List<Phylogeny> getPhylogenies() {
        return _phylogenies;
    }

    private Phylogeny[] getPhylogeniesAsArray() {
        final Phylogeny[] p = new Phylogeny[ getPhylogenies().size() ];
        for( int i = 0; i < getPhylogenies().size(); ++i ) {
            p[ i ] = getPhylogenies().get( i );
        }
        return p;
    }

    private List<String> getTaxlabels() {
        return _taxlabels;
    }

    private Map<String, String> getTranslateMap() {
        return _translate_map;
    }

    private boolean isIgnoreQuotes() {
        return _ignore_quotes_in_nh_data;
    }

    private boolean isReplaceUnderscores() {
        return _replace_underscores;
    }

    public Phylogeny[] parse() throws IOException, NHXFormatException {
        reset();
        final BufferedReader reader = ParserUtils.createReader( getNexusSource() );
        String line;
        String name = "";
        StringBuffer nhx = new StringBuffer();
        final StringBuffer translate_sb = new StringBuffer();
        boolean in_trees_block = false;
        boolean in_taxalabels = false;
        boolean in_translate = false;
        boolean in_tree = false;
        boolean rooted_info_present = false;
        boolean is_rooted = false;
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( ( line.length() > 0 ) && !line.startsWith( "#" ) && !line.startsWith( ">" ) ) {
                line = ForesterUtil.collapseWhiteSpace( line );
                line = removeWhiteSpaceBeforeSemicolon( line );
                final String line_lc = line.toLowerCase();
                if ( line_lc.startsWith( begin_trees ) ) {
                    in_trees_block = true;
                    in_taxalabels = false;
                    in_translate = false;
                }
                else if ( line_lc.startsWith( taxlabels ) ) {
                    in_trees_block = false;
                    in_taxalabels = true;
                    in_translate = false;
                }
                else if ( line_lc.startsWith( translate ) ) {
                    in_taxalabels = false;
                    in_translate = true;
                }
                else if ( in_trees_block ) {
                    if ( line_lc.startsWith( end ) || line_lc.startsWith( endblock ) ) {
                        in_trees_block = false;
                        in_tree = false;
                        in_translate = false;
                        if ( nhx.length() > 0 ) {
                            createPhylogeny( name, nhx, rooted_info_present, is_rooted );
                            nhx = new StringBuffer();
                            name = "";
                            rooted_info_present = false;
                            is_rooted = false;
                        }
                    }
                    else if ( line_lc.startsWith( tree ) || ( line_lc.startsWith( utree ) ) ) {
                        if ( nhx.length() > 0 ) {
                            createPhylogeny( name, nhx, rooted_info_present, is_rooted );
                            nhx = new StringBuffer();
                            name = "";
                            rooted_info_present = false;
                            is_rooted = false;
                        }
                        in_tree = true;
                        nhx.append( line.substring( line.indexOf( '=' ) ) );
                        final Matcher name_matcher = TREE_NAME_PATTERN.matcher( line );
                        if ( name_matcher.matches() ) {
                            name = name_matcher.group( 1 );
                            name = name.replaceAll( "['\"]+", "" );
                        }
                        final Matcher rootedness_matcher = ROOTEDNESS_PATTERN.matcher( line );
                        if ( rootedness_matcher.matches() ) {
                            final String s = rootedness_matcher.group( 1 );
                            line = line.replaceAll( "\\[\\&.\\]", "" );
                            rooted_info_present = true;
                            if ( s.toUpperCase().equals( "R" ) ) {
                                is_rooted = true;
                            }
                        }
                    }
                    else if ( in_tree && !in_translate ) {
                        nhx.append( line );
                    }
                    if ( !in_translate && !line_lc.startsWith( end ) && !line_lc.startsWith( endblock )
                            && line_lc.endsWith( ";" ) ) {
                        in_tree = false;
                        in_translate = false;
                        createPhylogeny( name, nhx, rooted_info_present, is_rooted );
                        nhx = new StringBuffer();
                        name = "";
                        rooted_info_present = false;
                        is_rooted = false;
                    }
                }
                if ( in_taxalabels ) {
                    if ( line_lc.startsWith( end ) || line_lc.startsWith( endblock ) ) {
                        in_taxalabels = false;
                    }
                    else {
                        final String[] labels = line.split( "\\s+" );
                        for( String label : labels ) {
                            if ( !label.toLowerCase().equals( taxlabels ) ) {
                                if ( label.endsWith( ";" ) ) {
                                    in_taxalabels = false;
                                    label = label.substring( 0, label.length() - 1 );
                                }
                                if ( label.length() > 0 ) {
                                    getTaxlabels().add( label );
                                }
                            }
                        }
                    }
                }
                if ( in_translate ) {
                    if ( line_lc.startsWith( end ) || line_lc.startsWith( endblock ) ) {
                        in_translate = false;
                    }
                    else {
                        translate_sb.append( " " );
                        translate_sb.append( line.trim() );
                        if ( line.endsWith( ";" ) ) {
                            in_translate = false;
                            setTranslateKeyValuePairs( translate_sb );
                        }
                    }
                }
            }
        }
        if ( nhx.length() > 0 ) {
            createPhylogeny( name, nhx, rooted_info_present, is_rooted );
        }
        return getPhylogeniesAsArray();
    }

    private void reset() {
        setPhylogenies( new ArrayList<Phylogeny>() );
        setTaxlabels( new ArrayList<String>() );
        setTranslateMap( new HashMap<String, String>() );
    }

    public void setIgnoreQuotes( final boolean ignore_quotes_in_nh_data ) {
        _ignore_quotes_in_nh_data = ignore_quotes_in_nh_data;
    }

    private void setPhylogenies( final ArrayList<Phylogeny> phylogenies ) {
        _phylogenies = phylogenies;
    }

    public void setReplaceUnderscores( final boolean replace_underscores ) {
        _replace_underscores = replace_underscores;
    }

    public void setSource( final Object nexus_source ) throws PhylogenyParserException, IOException {
        if ( nexus_source == null ) {
            throw new PhylogenyParserException( getClass() + ": attempt to parse null object." );
        }
        _nexus_source = nexus_source;
    }

    private void setTaxlabels( final List<String> taxlabels ) {
        _taxlabels = taxlabels;
    }

    private void setTranslateKeyValuePairs( final StringBuffer translate_sb ) throws IOException {
        String s = translate_sb.toString().trim();
        if ( s.endsWith( ";" ) ) {
            s = s.substring( 0, s.length() - 1 ).trim();
        }
        for( final String pair : s.split( "," ) ) {
            final String[] kv = pair.trim().split( "\\s+" );
            if ( ( kv.length < 2 ) || ( kv.length > 3 ) ) {
                throw new IOException( "ill formatted translate values: " + translate_sb );
            }
            if ( ( kv.length == 3 ) && !kv[ 0 ].toLowerCase().trim().equals( translate ) ) {
                throw new IOException( "ill formatted translate values: " + translate_sb );
            }
            String key = "";
            String value = "";
            if ( kv.length == 3 ) {
                key = kv[ 1 ];
                value = kv[ 2 ];
            }
            else {
                key = kv[ 0 ];
                value = kv[ 1 ];
            }
            if ( value.endsWith( ";" ) ) {
                value = value.substring( 0, value.length() - 1 );
            }
            getTranslateMap().put( key, value );
        }
    }

    private void setTranslateMap( final Map<String, String> translate_map ) {
        _translate_map = translate_map;
    }

    private static String removeWhiteSpaceBeforeSemicolon( final String s ) {
        return s.replaceAll( "\\s+;", ";" );
    }
}
