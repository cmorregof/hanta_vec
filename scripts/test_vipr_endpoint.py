#!/usr/bin/env python3
"""
Test ViPR web download endpoint directly as requested.
Report explicitly if unreachable rather than falling back silently.
"""

import requests
import time

def test_vipr_endpoint():
    """Test ViPR web download endpoint explicitly."""

    print("="*70)
    print("TESTING ViPR WEB DOWNLOAD ENDPOINT")
    print("="*70)

    endpoint = "https://www.viprbrc.org/brc/vipr_genome_search.spg"

    print(f"Testing endpoint: {endpoint}")

    try:
        # Test basic connectivity
        print("1. Testing basic connectivity...")
        response = requests.get(endpoint, timeout=30)

        print(f"   HTTP Status: {response.status_code}")
        print(f"   Response size: {len(response.text)} characters")

        if response.status_code == 200:
            print("   ✓ Endpoint is reachable")

            # Check if it's the correct ViPR page
            if "vipr" in response.text.lower() or "virus pathogen" in response.text.lower():
                print("   ✓ Appears to be ViPR website")
            else:
                print("   ⚠️ May not be the correct ViPR page")

        else:
            print(f"   ❌ HTTP error: {response.status_code}")
            return False

    except requests.exceptions.Timeout:
        print("   ❌ TIMEOUT: Endpoint is unreachable (30s timeout)")
        return False

    except requests.exceptions.ConnectionError:
        print("   ❌ CONNECTION ERROR: Cannot reach endpoint")
        return False

    except Exception as e:
        print(f"   ❌ ERROR: {e}")
        return False

    # Test form submission
    print("\n2. Testing form submission...")

    try:
        form_data = {
            'method': 'SubmitForm',
            'family': 'Hantaviridae',
            'download_format': 'fasta'
        }

        response = requests.post(endpoint, data=form_data, timeout=60)

        print(f"   Form submission status: {response.status_code}")

        if response.status_code == 200:
            print("   ✓ Form submission successful")

            # Check if response contains sequence data
            if ">" in response.text and ("hanta" in response.text.lower() or "virus" in response.text.lower()):
                print("   ✓ Response appears to contain FASTA sequences")

                # Count potential sequences
                sequence_count = response.text.count(">")
                print(f"   Potential sequences found: {sequence_count}")

                if sequence_count > 0:
                    print("   ✅ ViPR endpoint is FUNCTIONAL for Hantaviridae download")
                    return True
                else:
                    print("   ❌ No sequences found in response")
                    return False
            else:
                print("   ❌ Response does not contain expected FASTA format")
                print(f"   Response preview: {response.text[:200]}...")
                return False
        else:
            print(f"   ❌ Form submission failed: {response.status_code}")
            return False

    except requests.exceptions.Timeout:
        print("   ❌ TIMEOUT: Form submission timed out (60s)")
        return False

    except Exception as e:
        print(f"   ❌ ERROR during form submission: {e}")
        return False

def main():
    """Main function to test ViPR endpoint."""

    is_functional = test_vipr_endpoint()

    print(f"\n{'='*70}")
    print("ViPR ENDPOINT TEST RESULTS")
    print(f"{'='*70}")

    if is_functional:
        print("✅ ViPR web download endpoint is FUNCTIONAL")
        print("   - Endpoint is reachable")
        print("   - Form submission works")
        print("   - FASTA download available for Hantaviridae")
        print("\n💡 ViPR can be used as alternative to NCBI for sequence fetching")
    else:
        print("❌ ViPR web download endpoint is NOT FUNCTIONAL")
        print("   - Either unreachable or not working as expected")
        print("   - Continue using NCBI as primary source")
        print("\n⚠️  EXPLICIT REPORT: ViPR endpoint cannot be used from this environment")

    return 0 if is_functional else 1

if __name__ == "__main__":
    exit(main())