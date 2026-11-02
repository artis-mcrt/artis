// Fill the version selector of the sidebar with the versions that Read the Docs holds.
// The reader then goes from the documentation of the main branch to the documentation of a release.
// https://docs.readthedocs.com/platform/stable/addons.html describes the event and the data.
document.addEventListener("readthedocs-addons-data-ready", (event) => {
  const selector = document.querySelector(".version-selector");
  if (selector === null) {
    return;
  }
  const versions = event.detail.data().versions;
  const list = selector.querySelector("select");
  for (const version of versions.active) {
    const entry = document.createElement("option");
    entry.value = version.urls.documentation;
    entry.text = version.slug;
    entry.selected = version.slug === versions.current.slug;
    list.add(entry);
  }
  list.addEventListener("change", () => {
    window.location.href = list.value;
  });
  selector.hidden = false;
});
