export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
export HOME=/user/$USER

voms-proxy-init --voms cms -valid 192:00
