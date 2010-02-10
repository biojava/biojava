// $Id: PhyloXmlParser.java,v 1.11 2009/11/23 23:53:23 cmzmasek Exp $
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

package org.forester.io.parsers.phyloxml;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
import java.net.URL;
import java.util.Date;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.PhylogenyParserException;
import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.SAXNotSupportedException;
import org.xml.sax.SAXParseException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;

public class PhyloXmlParser implements PhylogenyParser {

    final public static String   JAXP_SCHEMA_LANGUAGE                       = "http://java.sun.com/xml/jaxp/properties/schemaLanguage";
    final public static String   W3C_XML_SCHEMA                             = "http://www.w3.org/2001/XMLSchema";
    final public static String   JAXP_SCHEMA_SOURCE                         = "http://java.sun.com/xml/jaxp/properties/schemaSource";
    final public static String   SAX_FEATURES_VALIDATION                    = "http://xml.org/sax/features/validation";
    final public static String   APACHE_FEATURES_VALIDATION_SCHEMA          = "http://apache.org/xml/features/validation/schema";
    final public static String   APACHE_FEATURES_VALIDATION_SCHEMA_FULL     = "http://apache.org/xml/features/validation/schema-full-checking";
    final public static String   APACHE_PROPERTIES_SCHEMA_EXTERNAL_LOCATION = "http://apache.org/xml/properties/schema/external-schemaLocation";
    final static private boolean TIME                                       = false;
    private Object               _source;
    private boolean              _valid;
    private boolean              _zipped_inputstream;
    private int                  _error_count;
    private int                  _warning_count;
    private String               _schema_location;
    private StringBuffer         _error_messages;
    private StringBuffer         _warning_messages;

    public static PhyloXmlParser createPhyloXmlParserXsdValidating() {
        final PhyloXmlParser xml_parser = new PhyloXmlParser();
        final ClassLoader cl = PhyloXmlParser.class.getClassLoader();
        final URL xsd_url = cl.getResource( ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE );
        if ( xsd_url != null ) {
            xml_parser.setValidateAgainstSchema( xsd_url.toString() );
        }
        else {
            throw new IllegalStateException( "failed to get URL for phyloXML XSD from jar file from ["
                    + ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE + "]" );
        }
        return xml_parser;
    }

    public PhyloXmlParser() {
        init();
        reset();
    }

    public int getErrorCount() {
        return _error_count;
    }

    public StringBuffer getErrorMessages() {
        return _error_messages;
    }

    private Reader getReaderFromZipFile() throws IOException {
        Reader reader = null;
        final ZipFile zip_file = new ZipFile( getSource().toString() );
        final Enumeration<?> zip_file_entries = zip_file.entries();
        while ( zip_file_entries.hasMoreElements() ) {
            final ZipEntry zip_file_entry = ( ZipEntry ) zip_file_entries.nextElement();
            if ( !zip_file_entry.isDirectory() && ( zip_file_entry.getSize() > 0 ) ) {
                final InputStream is = zip_file.getInputStream( zip_file_entry );
                reader = new InputStreamReader( is );
                break;
            }
        }
        return reader;
    }

    private String getSchemaLocation() {
        return _schema_location;
    }

    private Object getSource() {
        return _source;
    }

    public int getWarningCount() {
        return _warning_count;
    }

    public StringBuffer getWarningMessages() {
        return _warning_messages;
    }

    private void init() {
        setZippedInputstream( false );
    }

    public boolean isValid() {
        return _valid;
    }

    private boolean isZippedInputstream() {
        return _zipped_inputstream;
    }

    public Phylogeny[] parse() throws IOException, PhylogenyParserException {
        reset();
        final PhyloXmlHandler handler = new PhyloXmlHandler();
        final SAXParserFactory factory = SAXParserFactory.newInstance();
        factory.setNamespaceAware( true );
        try {
            if ( !ForesterUtil.isEmpty( getSchemaLocation() ) ) {
                factory.setFeature( SAX_FEATURES_VALIDATION, true );
                factory.setFeature( APACHE_FEATURES_VALIDATION_SCHEMA, true );
                factory.setFeature( APACHE_FEATURES_VALIDATION_SCHEMA_FULL, true );
            }
        }
        catch ( final SAXNotRecognizedException e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "sax not recognized exception: " + e.getLocalizedMessage() );
        }
        catch ( final SAXNotSupportedException e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "sax not supported exception: " + e.getLocalizedMessage() );
        }
        catch ( final ParserConfigurationException e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "parser configuration exception: " + e.getLocalizedMessage() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "error while configuring sax parser: " + e.getLocalizedMessage() );
        }
        try {
            final SAXParser parser = factory.newSAXParser();
            if ( !ForesterUtil.isEmpty( getSchemaLocation() ) ) {
                parser.setProperty( JAXP_SCHEMA_LANGUAGE, W3C_XML_SCHEMA );
                parser.setProperty( JAXP_SCHEMA_SOURCE, getSchemaLocation() );
                parser.setProperty( APACHE_PROPERTIES_SCHEMA_EXTERNAL_LOCATION, getSchemaLocation() );
            }
            final XMLReader xml_reader = parser.getXMLReader();
            xml_reader.setContentHandler( handler );
            xml_reader.setErrorHandler( new PhyloXmlParserErrorHandler() );
            long start_time = 0;
            if ( TIME ) {
                start_time = new Date().getTime();
            }
            if ( getSource() instanceof File ) {
                if ( !getSource().toString().toLowerCase().endsWith( ".zip" ) ) {
                    xml_reader.parse( new InputSource( new FileReader( ( File ) getSource() ) ) );
                }
                else {
                    final Reader reader = getReaderFromZipFile();
                    if ( reader == null ) {
                        throw new PhylogenyParserException( "zip file \"" + getSource()
                                + "\" appears not to contain any entries" );
                    }
                    xml_reader.parse( new InputSource( reader ) );
                }
            }
            else if ( getSource() instanceof InputSource ) {
                xml_reader.parse( ( InputSource ) getSource() );
            }
            else if ( getSource() instanceof InputStream ) {
                if ( !isZippedInputstream() ) {
                    final InputStream is = ( InputStream ) getSource();
                    final Reader reader = new InputStreamReader( is );
                    xml_reader.parse( new InputSource( reader ) );
                }
                else {
                    final ZipInputStream zip_is = new ZipInputStream( ( InputStream ) getSource() );
                    zip_is.getNextEntry();
                    final Reader reader = new InputStreamReader( zip_is );
                    if ( reader == null ) {
                        throw new PhylogenyParserException( "zip input stream \"" + getSource()
                                + "\" appears not to contain any (phyloXML) data" );
                    }
                    xml_reader.parse( new InputSource( reader ) );
                }
            }
            else if ( getSource() instanceof String ) {
                final File file = new File( getSource().toString() );
                final Reader reader = new FileReader( file );
                xml_reader.parse( new InputSource( reader ) );
            }
            else if ( getSource() instanceof StringBuffer ) {
                final StringReader string_reader = new StringReader( getSource().toString() );
                xml_reader.parse( new InputSource( string_reader ) );
            }
            else {
                throw new PhylogenyParserException( "phyloXML parser: attempt to parse object of unsupported type: \""
                        + getSource().getClass() + "\"" );
            }
            if ( TIME ) {
                System.out.println( "[TIME] phyloXML parsing: " + ( new Date().getTime() - start_time ) + "ms." );
            }
        }
        catch ( final SAXException sax_exception ) {
            throw new PhylogenyParserException( "failed to parse [" + getSource() + "]: "
                    + sax_exception.getLocalizedMessage() );
        }
        catch ( final ParserConfigurationException parser_config_exception ) {
            throw new PhylogenyParserException( "failed to parse [" + getSource()
                    + "]. Problem with XML parser configuration: " + parser_config_exception.getLocalizedMessage() );
        }
        catch ( final IOException e ) {
            throw new PhylogenyParserException( "problem with input source: " + e.getLocalizedMessage() );
        }
        catch ( final Exception e ) {
            throw new PhylogenyParserException( e.getLocalizedMessage() );
        }
        catch ( final Error err ) {
            err.printStackTrace();
            throw new PhylogenyParserException( "severe error: " + err.getLocalizedMessage() );
        }
        final Phylogeny[] ps = new Phylogeny[ handler.getPhylogenies().size() ];
        int i = 0;
        for( final Phylogeny phylogeny : handler.getPhylogenies() ) {
            ps[ i++ ] = phylogeny;
        }
        return ps;
    }

    private void reset() {
        _valid = true;
        _error_count = 0;
        _warning_count = 0;
        _error_messages = new StringBuffer();
        _warning_messages = new StringBuffer();
    }

    public void setSource( final Object source ) {
        _source = source;
    }

    public void setValidateAgainstSchema( final String schema_location ) {
        _schema_location = schema_location;
    }

    public void setZippedInputstream( final boolean zipped_inputstream ) {
        _zipped_inputstream = zipped_inputstream;
    }

    private class PhyloXmlParserErrorHandler extends DefaultHandler {

        @Override
        public void error( final SAXParseException e ) {
            ++_error_count;
            _valid = false;
            throw new PhyloXmlException( "phyloXML error at line " + e.getLineNumber() + ": \n"
                    + e.getLocalizedMessage() );
        }

        @Override
        public void fatalError( final SAXParseException e ) {
            ++_error_count;
            _valid = false;
            throw new PhyloXmlException( "fatal XML error at line " + e.getLineNumber() + ": \n"
                    + e.getLocalizedMessage() );
        }

        @Override
        public void warning( final SAXParseException e ) {
            ++_warning_count;
            if ( _error_messages.length() > 1 ) {
                _error_messages.append( ForesterUtil.LINE_SEPARATOR );
            }
            _warning_messages.append( "[line: " + e.getLineNumber() + "] " + e.getMessage() );
        }
    }
}
