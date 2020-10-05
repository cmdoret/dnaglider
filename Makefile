
GOFILES = $(shell find . -name '*.go')
VERSION="0.0.1"

default: release 

build/:
	mkdir -p $@

.PHONY: deps
deps:
	go get -d github.com/cmdoret/dnaglider/dnaglider

build/dnaglider-windows.exe: $(GOFILES) build/ deps
	GOOS=windows GOARCH=amd64 CGO_ENABLED=0 go build -ldflags "-X main.VERSION=$(VERSION)" -o $@ ./dnaglider/

build/dnaglider-osx: $(GOFILES) build/ deps
	GOOS=darwin GOARCH=amd64 CGO_ENABLED=0 go build -ldflags "-X main.VERSION=$(VERSION)" -o $@ ./dnaglider/

build/dnaglider-linux: $(GOFILES) build/ deps
	GOOS=linux GOARCH=amd64 CGO_ENABLED=0 go build -ldflags "-X main.VERSION=$(VERSION)" -o $@ ./dnaglider/

.PHONY: release
release: build/dnaglider-linux build/dnaglider-windows.exe build/dnaglider-osx

.PHONY: test
test: deps
	cd dnaglider/pkg ; go test
