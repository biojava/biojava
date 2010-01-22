// $Id: XmlElement.java,v 1.8 2009/10/30 02:53:54 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.HashMap;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.util.ForesterUtil;
import org.xml.sax.Attributes;

public class XmlElement {

    public final static boolean           DEBUG = false;
    private final String                  _namespaceUri;
    private final String                  _localName;
    private final String                  _qualifiedName;
    private String                        _value;
    private final HashMap<String, String> _attributes;
    private final ArrayList<XmlElement>   _childElements;
    private XmlElement                    _parent;

    public XmlElement( final String namespaceUri,
                       final String localName,
                       final String qualifiedName,
                       final Attributes attributes ) {
        _namespaceUri = namespaceUri;
        _localName = localName;
        _qualifiedName = qualifiedName;
        if ( attributes != null ) {
            _attributes = new HashMap<String, String>( attributes.getLength() );
            for( int i = 0; i < attributes.getLength(); ++i ) {
                getAttributes().put( new String( attributes.getQName( i ) ), new String( attributes.getValue( i ) ) );
            }
        }
        else {
            _attributes = new HashMap<String, String>();
        }
        _childElements = new ArrayList<XmlElement>();
        _parent = null;
    }

    public void addChildElement( final XmlElement element ) {
        element.setParent( this );
        getChildElements().add( element );
    }

    public void appendValue( final String value ) {
        _value = _value + value;
    }

    public String getAttribute( final String attribute_name ) {
        if ( !isHasAttribute( attribute_name ) ) {
            throw new IllegalArgumentException( "no attribute named [" + attribute_name + "] present in element ["
                    + getQualifiedName() + "]" );
        }
        return getAttributes().get( attribute_name );
    }

    public HashMap<String, String> getAttributes() {
        return _attributes;
    }

    public XmlElement getChildElement( final int i ) {
        if ( ( i < 0 ) || ( i >= getNumberOfChildElements() ) ) {
            throw new IllegalArgumentException( "attempt to get child element with index " + i + " for element with "
                    + getNumberOfChildElements() + " child elements" );
        }
        return getChildElements().get( i );
    }

    ArrayList<XmlElement> getChildElements() {
        return _childElements;
    }

    String getLocalName() {
        return _localName;
    }

    String getNamespaceUri() {
        return _namespaceUri;
    }

    public int getNumberOfChildElements() {
        return getChildElements().size();
    }

    public XmlElement getParent() {
        return _parent;
    }

    public String getQualifiedName() {
        return _qualifiedName;
    }

    XmlElement getRoot() {
        XmlElement e = this;
        while ( e.getParent() != null ) {
            e = e.getParent();
        }
        return e;
    }

    public boolean getValueAsBoolean() throws PhylogenyParserException {
        boolean b = false;
        try {
            b = ( new Boolean( getValueAsString() ) ).booleanValue();
        }
        catch ( final NumberFormatException ex ) {
            throw new PhylogenyParserException( "attempt to parse [" + getValueAsString() + "] into boolean, in "
                    + toString() );
        }
        return b;
    }

    public double getValueAsDouble() throws PhylogenyParserException {
        double d = 0.0;
        try {
            d = Double.parseDouble( getValueAsString() );
        }
        catch ( final NumberFormatException ex ) {
            throw new PhylogenyParserException( "attempt to parse [" + getValueAsString() + "] into double, in "
                    + toString() );
        }
        return d;
    }

    public int getValueAsInt() throws PhylogenyParserException {
        int i = 0;
        try {
            i = Integer.parseInt( getValueAsString() );
        }
        catch ( final NumberFormatException ex ) {
            throw new PhylogenyParserException( "attempt to parse [" + getValueAsString() + "] into integer, in "
                    + toString() );
        }
        return i;
    }

    public String getValueAsString() {
        if ( _value == null ) {
            return "";
        }
        return _value.replaceAll( "\\s+", " " ).trim();
    }

    public boolean isHasAttribute( final String attribute_name ) {
        return getAttributes().containsKey( attribute_name );
    }

    public boolean isHasValue() {
        return !ForesterUtil.isEmpty( _value );
    }

    void setParent( final XmlElement parent ) {
        _parent = parent;
    }

    /**
     * [Careful, this does not call "new String(...)"]
     * 
     * @param value
     */
    public void setValue( final String value ) {
        _value = value;
        if ( XmlElement.DEBUG ) {
            System.out.println();
            System.out.println( "Value is \"" + value + "\" for" );
            System.out.println( "Local name     = " + getLocalName() );
            System.out.println( "Qualified name = " + getQualifiedName() );
            System.out.println( "Namespace URI  = " + getNamespaceUri() );
            System.out.print( "Attributes     : " );
            for( final String string : getAttributes().keySet() ) {
                final String key = string;
                System.out.print( key + " = \"" + getAttributes().get( key ) + "\"  " );
            }
            System.out.println();
            System.out.println();
        }
    }

    @Override
    public String toString() {
        if ( getParent() != null ) {
            return "\"" + getQualifiedName() + "\" [value: " + getValueAsString() + ", parent element: \""
                    + getParent().getQualifiedName() + "\"]";
        }
        return "\"" + getQualifiedName() + "\" [value: " + getValueAsString() + "]";
    }
}
