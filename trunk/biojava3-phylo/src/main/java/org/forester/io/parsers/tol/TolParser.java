// $Id: TolParser.java,v 1.10 2009/10/26 23:29:40 cmzmasek Exp $
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

package org.forester.io.parsers.tol;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;
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
import org.forester.util.ForesterUtil;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.SAXNotSupportedException;
import org.xml.sax.SAXParseException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;

public class TolParser implements PhylogenyParser {

    final public static String JAXP_SCHEMA_LANGUAGE                       = "http://java.sun.com/xml/jaxp/properties/schemaLanguage";
    final public static String W3C_XML_SCHEMA                             = "http://www.w3.org/2001/XMLSchema";
    final public static String JAXP_SCHEMA_SOURCE                         = "http://java.sun.com/xml/jaxp/properties/schemaSource";
    final public static String SAX_FEATURES_VALIDATION                    = "http://xml.org/sax/features/validation";
    final public static String APACHE_FEATURES_VALIDATION_SCHEMA          = "http://apache.org/xml/features/validation/schema";
    final public static String APACHE_FEATURES_VALIDATION_SCHEMA_FULL     = "http://apache.org/xml/features/validation/schema-full-checking";
    final public static String APACHE_PROPERTIES_SCHEMA_EXTERNAL_LOCATION = "http://apache.org/xml/properties/schema/external-schemaLocation";
    private Object             _source;
    private boolean            _valid;
    private boolean            _zipped_inputstream;
    private int                _error_count;
    private int                _warning_count;
    private String             _schema_location;
    private StringBuffer       _error_messages;
    private StringBuffer       _warning_messages;

    public TolParser() {
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
        final TolXmlHandler handler = new TolXmlHandler();
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
            throw new PhylogenyParserException( "sax not recognized exception: " + e.getMessage() );
        }
        catch ( final SAXNotSupportedException e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "sax not supported exception: " + e.getMessage() );
        }
        catch ( final ParserConfigurationException e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "parser _configuration exception: " + e.getMessage() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "error while configuring sax parser: " + e.getMessage() );
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
            xml_reader.setErrorHandler( new TolParserErrorHandler() );
            if ( getSource() instanceof File ) {
                if ( !getSource().toString().toLowerCase().endsWith( ".zip" ) ) {
                    xml_reader.parse( new InputSource( new FileReader( ( File ) getSource() ) ) );
                }
                else {
                    final Reader reader = getReaderFromZipFile();
                    if ( reader == null ) {
                        throw new PhylogenyParserException( "Zip file \"" + getSource()
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
                        throw new PhylogenyParserException( "Zip input stream \"" + getSource()
                                + "\" appears not to contain any data" );
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
                throw new PhylogenyParserException( "attempt to parse object of unsupported type: \""
                        + getSource().getClass() + "\"" );
            }
        }
        catch ( final SAXException sax_exception ) {
            throw new PhylogenyParserException( "Failed to parse [" + getSource() + "]: " + sax_exception.getMessage() );
        }
        catch ( final ParserConfigurationException parser_config_exception ) {
            throw new PhylogenyParserException( "Failed to parse [" + getSource()
                    + "] Problem with xml parser _configuration: " + parser_config_exception.getMessage() );
        }
        catch ( final IOException e ) {
            throw new PhylogenyParserException( "Problem with input source [" + getSource() + "]: \n" + e.getMessage() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            throw new PhylogenyParserException( "Failed to parse [" + getSource() + "]: " + e.getMessage() );
        }
        catch ( final Error err ) {
            err.printStackTrace();
            throw new PhylogenyParserException( "Severe error: " + err.getMessage() );
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

    private class TolParserErrorHandler extends DefaultHandler {

        @Override
        public void error( final SAXParseException e ) {
            ++_error_count;
            _valid = false;
            throw new RuntimeException( "XML error at line " + e.getLineNumber() + ": \n" + e.getMessage() );
        }

        @Override
        public void fatalError( final SAXParseException e ) {
            ++_error_count;
            _valid = false;
            throw new RuntimeException( "Fatal XML error at line " + e.getLineNumber() + ": \n" + e.getMessage() );
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