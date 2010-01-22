// $Id: Identifier.java,v 1.41 2009/11/10 21:16:37 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class Identifier implements PhylogenyData {

    final String _value;
    final String _provider;

    public Identifier() {
        _value = "";
        _provider = "";
    }

    public Identifier( final String value ) {
        _value = value;
        _provider = "";
    }

    public Identifier( final String value, final String provider ) {
        _value = value;
        _provider = provider;
    }

    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue() );
    }

    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getProvider() ) ) {
            sb.append( "[" );
            sb.append( getProvider() );
            sb.append( "] " );
        }
        sb.append( getValue() );
        return sb;
    }

    public PhylogenyData copy() {
        return new Identifier( new String( getValue() ), new String( getProvider() ) );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            return false;
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return isEqual( ( Identifier ) o );
        }
    }

    public String getProvider() {
        return _provider;
    }

    public String getValue() {
        return _value;
    }

    @Override
    public int hashCode() {
        if ( getProvider() != null ) {
            return ( getProvider() + getValue() ).hashCode();
        }
        return getValue().hashCode();
    }

    public boolean isEqual( final PhylogenyData data ) {
        if ( this == data ) {
            return true;
        }
        if ( ( data == null ) || ( getValue() == null ) ) {
            return false;
        }
        final Identifier a = ( Identifier ) data;
        if ( ( getProvider() != null ) && ( a.getProvider() != null ) ) {
            return ( a.getValue().equals( getValue() ) && a.getProvider().equals( getProvider() ) );
        }
        return ( a.getValue().equals( getValue() ) );
    }

    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( ":" );
        sb.append( NHXtags.NODE_IDENTIFIER );
        sb.append( ForesterUtil.replaceIllegalNhxCharacters( getValue() ) );
        return sb;
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( !org.forester.util.ForesterUtil.isEmpty( getProvider() ) ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.IDENTIFIER,
                                             getValue(),
                                             PhyloXmlMapping.IDENTIFIER_PROVIDER_ATTR,
                                             getProvider(),
                                             indentation );
        }
        else {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.IDENTIFIER, getValue(), indentation );
        }
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
